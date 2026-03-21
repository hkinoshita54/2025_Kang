# cellchat of primary tumor sample
# only 3 clusters: Tumor, Fibro, EC
# use updated cellchatdb from 2025_Kumar_2
# 2026-03-10

# Make directories ----
analysis_step <- "300_cellchat"
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

## load libraries
library(CellChat)
library(patchwork)
library(Seurat)
source("Rscripts/theme_panel.R")
options(stringsAsFactors = FALSE)
library(forcats)

# create seurat obj for cellchat ----
## Load seurat objects, all cells and annotated cells
seu_all <- readRDS("RDSfiles/seu_020_all.RDS")

epi <- readRDS("RDSfiles/seu_030_epi.RDS")
tum <- subset(epi, subset = celltype2 == "Tumor")
tum_names <- Cells(tum)

fib <- subset(seu_all, subset = celltype1 == "Fibro" & tissue_type == "T")
fib_names <- Cells(fib)

ec <- subset(seu_all, subset = celltype1 == "EC" & tissue_type == "T")
ec_names <- Cells(ec)

seu <- subset(seu_all, cells = c(tum_names, fib_names, ec_names))

# check cell numbers ----
seu$pt_id <- droplevels(seu$pt_id)
seu$celltype1 <- droplevels(seu$celltype1)
table(seu$pt_id, seu$celltype1)

## remove samples where any of the cell types has < 20 cells
tab <- table(seu$pt_id, seu$celltype1)
keep <- rownames(tab)[apply(tab, 1, function(x) all(x >= 20))]
seu <- subset(seu, subset = pt_id %in% keep)
seu$pt_id <- droplevels(seu$pt_id)
table(seu$pt_id, seu$celltype1)

seu <- DietSeurat(seu)
seu_all <- seu

# CellChat analysis per patient ----
load("RDSfiles/CellChatDB.new.RData")
CellChatDB.use <- subsetDB(CellChatDB.new, search = c("Secreted Signaling"), key = "annotation")
seu_all <- NormalizeData(seu_all) # need normalized count for cellchat
pts <- levels(seu_all$pt_id)
cellchat_list <- list()

for (pid in pts){
  cat("Processing:", pid, "\n")
  seu <- subset(seu_all, subset = pt_id == pid)
  cellchat <- createCellChat(object = seu, group.by = "celltype1", assay = "RNA")
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cellchat_list[[pid]] <- cellchat
}

names(cellchat_list) <- pts

comm_list <- lapply(names(cellchat_list), function(pid) {
  df <- subsetCommunication(cellchat_list[[pid]])
  df$patient_ID <- pid
  df
})

comm_df <- do.call(rbind, comm_list)

comm_summary <- comm_df %>%
  group_by(source, target, ligand, receptor, pathway_name) %>%
  summarise(
    n_patients = n_distinct(patient_ID),
    mean_prob = mean(prob, na.rm = TRUE),
    median_prob = median(prob, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_patients), desc(mean_prob))

openxlsx2::write_xlsx(comm_df, file = file.path(res_path, "comm_df.xlsx"))
openxlsx2::write_xlsx(comm_summary, file = file.path(res_path, "comm_summary.xlsx"))

saveRDS(cellchat_list, file.path("RDSfiles", "cellchat_list_321.RDS"))

# visualization ----
df <- openxlsx2::read_xlsx("result/300_cellchat/comm_summary.xlsx")

df_plot <- df %>%
  filter(source == "Epi", target == "EC") %>%
  group_by(ligand) %>%
  arrange(desc(n_patients), desc(mean_prob), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    lr_label = paste0(ligand, " - ", receptor),
    pathway_name = factor(pathway_name)
  ) %>%
  arrange(desc(n_patients), desc(mean_prob)) %>%
  mutate(
    highlight = ifelse(ligand == "LRG1" & receptor == "ENG",
                       "LRG1-ENG", "Other")
  )
df_plot <- df_plot %>% slice_head(n=10)

ggplot(df_plot, aes(x = fct_reorder(lr_label, n_patients), y = n_patients, fill = highlight)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("Other" = "grey35", "LRG1-ENG" = "firebrick")) +
  scale_y_continuous(
    breaks = seq(0, max(df_plot$n_patients) + 1, by = 1),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = NULL,
    y = "n patients"
  ) +
  theme_panel() & NoLegend()
ggsave("Epi_EC.pdf", path = plot_path, width = 40, height = 30, units = "mm")
