library(tidyverse)
library(ggbeeswarm)
library(ggpmisc)

# custom x or y axis breaks
log10_minor_break = function (...) {
  function(x) {
    minx         = floor(min(log10(x), na.rm = T)) - 1
    maxx         = ceiling(max(log10(x), na.rm = T)) + 1
    n_major      = maxx - minx + 1
    major_breaks = seq(minx, maxx, by = 1)
    minor_breaks =
      rep(log10(seq(1, 9, by = 1)), times = n_major) +
      rep(major_breaks, each = 9)
    return(10 ^ (minor_breaks))
  }
}

# QC filtered values with NIST SRM 1950 values
quantities <-
  read_csv("LipidCreatorVal_QC-filtered_SAMPLES-NIST_uM_04022020v1.csv")
qtys <-
  quantities %>% tidyr::pivot_longer(
    cols = 5:length(colnames(.)),
    names_to = "lc_lipid",
    values_to = "lc_quantity"
  )
#group by lipid and subject and fix some of the names
#SubjectID,BLOCK,REPLICATE
qtys_by_subject <- qtys %>%
  rename(SubjectID = SampleID) %>%
  separate(REPLICATE,
           into = c("BLOCK", "REPLICATE"),
           sep = 1) %>%
  mutate(lc_lipid = gsub(" ?[abc]+$", "", lc_lipid)) %>%
  mutate(lc_lipid = gsub("[abc]+$", "", lc_lipid)) %>%
  mutate(lc_lipid = gsub(" ?\\[NL[0-9:-]+\\]", "", lc_lipid)) %>%
  mutate(lc_lipid = gsub("^TG ", "TAG ", lc_lipid)) %>%
  mutate(lc_lipid = gsub("^DG ", "DAG ", lc_lipid)) %>%
  mutate(lc_lipid = gsub("^Hex1Cer ", "HexCer ", lc_lipid)) %>%
  group_by(lc_lipid, SubjectID) %>%
  summarise("lc_mean" = mean(lc_quantity),
            "lc_stdev" = sd(lc_quantity)) %>%
  rename("Lipid" = "lc_lipid",
         "Consensus" = "lc_mean",
         "Uncertainty" = "lc_stdev") %>%
  mutate(size = 1) %>%
  ungroup() %>%
  mutate(Lipid = gsub("CE ", "ChE ", Lipid)) %>%
  mutate(Lipid = gsub("d18:0", "18:0;2", Lipid)) %>%
  mutate(Lipid = gsub("d18:1", "18:1;2", Lipid)) %>%
  mutate(Lipid = factor(Lipid)) %>%
  group_by(Lipid, SubjectID)

qtys_nist <- qtys_by_subject %>% filter(SubjectID == "NIST")
qtys_subjects <- qtys_by_subject %>% filter(SubjectID != "NIST")

# summarize quantities by lipid species
qtys_summary <- qtys %>% filter(SampleID != "NIST") %>%
  mutate(lc_lipid = gsub(" ?[abc]+$", "", lc_lipid)) %>%
  mutate(lc_lipid = gsub("[abc]+$", "", lc_lipid)) %>%
  mutate(lc_lipid = gsub(" ?\\[NL[0-9:-]+\\]", "", lc_lipid)) %>%
  mutate(lc_lipid = gsub("^TG ", "TAG ", lc_lipid)) %>%
  mutate(lc_lipid = gsub("^DG ", "DAG ", lc_lipid)) %>%
  mutate(lc_lipid = gsub("^Hex1Cer ", "HexCer ", lc_lipid)) %>%
  group_by(lc_lipid) %>%
  summarise(
    "lc_mean" = mean(lc_quantity),
    "lc_stdev" = sd(lc_quantity),
    "lc_std_error" = sd(lc_quantity) / sqrt(n())
  ) %>%
  mutate(lc_lipid = factor(lc_lipid))

qtys_nist_summary <- qtys %>% filter(SampleID == "NIST") %>%
  mutate(lc_lipid = gsub(" ?[abc]+$", "", lc_lipid)) %>%
  mutate(lc_lipid = gsub("[abc]+$", "", lc_lipid)) %>%
  mutate(lc_lipid = gsub(" ?\\[NL[0-9:-]+\\]", "", lc_lipid)) %>%
  mutate(lc_lipid = gsub("^TG ", "TAG ", lc_lipid)) %>%
  mutate(lc_lipid = gsub("^DG ", "DAG ", lc_lipid)) %>%
  mutate(lc_lipid = gsub("^Hex1Cer ", "HexCer ", lc_lipid)) %>%
  group_by(lc_lipid) %>%
  summarise(
    "lc_mean" = mean(lc_quantity),
    "lc_stdev" = sd(lc_quantity),
    "lc_std_error" = sd(lc_quantity) / sqrt(n())
  ) %>%
  mutate(lc_lipid = factor(lc_lipid))

# renaming and factor creation
qtys_summary_renamed <-
  qtys_summary %>%
  select(-lc_std_error) %>%
  rename("Consensus" = "lc_mean",
         "Uncertainty" = "lc_stdev",
         "Lipid" = "lc_lipid") %>%
  mutate(origin = "Mean Individuals") %>%
  mutate(Lipid = gsub("CE ", "ChE ", Lipid)) %>%
  mutate(Lipid = gsub("d18:0", "18:0;2", Lipid)) %>%
  mutate(Lipid = gsub("d18:1", "18:1;2", Lipid)) %>%
  mutate(Lipid = factor(Lipid))

qtys_nist_renamed <-
  qtys_nist_summary %>%
  select(-lc_std_error) %>%
  rename("Consensus" = "lc_mean",
         "Uncertainty" = "lc_stdev",
         "Lipid" = "lc_lipid") %>%
  mutate(origin = "NIST SRM 1950") %>%
  mutate(Lipid = gsub("CE ", "ChE ", Lipid)) %>%
  mutate(Lipid = gsub("d18:0", "18:0;2", Lipid)) %>%
  mutate(Lipid = gsub("d18:1", "18:1;2", Lipid)) %>%
  mutate(Lipid = factor(Lipid))

# names of target lipids
lipidNames <- data.frame(Lipid = qtys_summary_renamed$Lipid)

# lipids to show in the panel plot
panelLipids <- c(
  "ChE 16:1",
  "ChE 18:0",
  "ChE 18:1",
  "ChE 18:2",
  "ChE 18:3",
  "ChE 20:3",
  "ChE 22:6",
  "Cer 18:1;2/16:0",
  "Cer 18:1;2/18:0",
  "Cer 18:1;2/24:1",
  "HexCer 18:1;2/16:0",
  "HexCer 18:1;2/18:0",
  "HexCer 18:1;2/22:0",
  "HexCer 18:1;2/24:1",
  "LPC 16:0",
  "LPC 18:0",
  "LPC 18:1",
  "LPC 18:2",
  "LPC 20:4",
  "LPC 22:6",
  "LPC O-18:0a",
  "LPE 16:0",
  "LPE 18:0",
  "PC 34:1",
  "PC 34:2",
  "PC 36:1",
  "PC 36:2",
  "PC 36:3",
  "PC 36:4",
  "PC 36:5",
  "PC 38:6",
  "PC 40:6",
  "PE 34:2",
  "PE 36:2",
  "PE 36:3",
  "PE 36:4",
  "PE 38:4",
  "PE 38:5",
  "PE 40:6",
  "PI 36:2",
  "PI 36:3",
  "PI 36:4",
  "PI 38:4",
  "SM 34:1",
  "SM 34:2",
  "SM 36:1",
  "SM 36:2",
  "SM 38:1",
  "TAG 54:5",
  "TAG 56:6"
)

# create a common lipids data frame from Bowden et al. data and our measurements.
commonLipids <-
  bind_rows(qtys_summary_renamed, qtys_nist_renamed) %>%
  left_join(lipidNames, by = c("Lipid" = "Lipid")) %>%
  mutate(
    Lipid = as.factor(Lipid),
    origin = as.factor(origin),
    ymin = Consensus - Uncertainty,
    ymax = Consensus + Uncertainty
  ) %>%
  filter(Lipid %in% panelLipids) %>%
  mutate(Lipid = factor(Lipid, levels = panelLipids, ordered = TRUE))

# common lipid species with refactoring to have defined plot order
commonLipidSpecies <-
  qtys_subjects %>% filter(Lipid %in% panelLipids) %>% ungroup() %>%
  mutate(Lipid = factor(Lipid, levels = panelLipids, ordered = TRUE)) %>% group_by(Lipid)

#write_csv(x = commonLipids, path = "plasma-vs-NIST-common-lipids.csv")
cscalefill <- scale_fill_manual(values = c("#d95f02", "#7570b3"))
cscalestroke <- scale_color_manual(values = c("#d95f02", "#7570b3"))
# plot the facetted comparison plot, colored areas represent Uncertainty regions (std uncertainty for Bowden et al., std deviation for our measurements)
facetComparisonPlot <- ggplot(data = commonLipids) +
  geom_quasirandom(
    pch = 16,
    groupOnX = TRUE,
    method = "quasirandom",
    varwidth = FALSE,
    width = 0.3,
    bandwidth = 0.1,
    aes(x = Lipid,
        y = Consensus,
        size = as.factor(size)),
    #alpha = 1.0,
    fill = "black",
    color = "black",
    data = commonLipidSpecies
  ) +
  geom_hline(aes(
    yintercept = Consensus,
    color = origin,
    linetype = origin
  ),
  size = 1.1) +
  geom_rect(
    aes(ymin = ymin, ymax = ymax, fill = origin),
    xmin = -Inf,
    xmax = Inf,
    alpha = 0.1,
    show.legend = FALSE
  ) +
  expand_limits(y = 0) +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_rect(
      fill = "white",
      size = 0.5,
      linetype = "solid",
      colour = "white",
      
    ),
    text = element_text(size=12),
    strip.text = element_text(size=8, face = "bold"),
    plot.caption = element_text(hjust = 1),
    axis.title.y = element_text(face="bold")
  ) + guides(
    colour = guide_legend(order = 1),
    linetype = guide_legend(order = 1),
    # alpha = guide_legend(override.aes = list(size = 2), order = 2),
    size = guide_legend(label = TRUE, override.aes = list(size = 2), order = 2)
  ) +
  facet_wrap(. ~ Lipid, scales = "free", ncol = 7) +
  cscalefill +
  cscalestroke +
  scale_size_manual(labels = c("Individuals"),
                    values = c(2)) +
  scale_linetype_manual(values=c("dotdash", "dashed"))+
  ylab(bquote("Concentration" ~ group("(",paste(mu, mol, ~L^-1), ")"))) +
  labs(caption = bquote("Shaded areas: " ~ mu %+-% sigma))
ggsave(
  "plasmaComparisonPlot.pdf",
  facetComparisonPlot,
  width = 11.69,
  height = 8.27,
  units = "in",
  useDingbats = FALSE
)
ggsave(
  "plasmaComparisonPlot.svg",
  facetComparisonPlot,
  width = 11.69,
  height = 8.27,
  units = "in"
)

# plot the concentration summary by lipid class overview
concentrationData <- qtys_subjects %>% mutate(Class = Lipid) %>%
  mutate(Class = gsub("ChE.*", "ChE", Class)) %>%
  mutate(Class = gsub("^Cer.*", "Cer", Class)) %>%
  mutate(Class = gsub("^HexCer.*", "HexCer", Class)) %>%
  mutate(Class = gsub("^Hex2Cer.*", "Hex2Cer", Class)) %>%
  mutate(Class = gsub("^Hex3Cer.*", "Hex3Cer", Class)) %>%
  mutate(Class = gsub("^LPC O-.*", "LPC-O", Class)) %>%
  mutate(Class = gsub("^LPE .*", "LPE", Class)) %>%
  mutate(Class = gsub("^LPC .*", "LPC", Class)) %>%
  mutate(Class = gsub("^LPI .*", "LPI", Class)) %>%
  mutate(Class = gsub("^PC O-.*", "PC-O", Class)) %>%
  mutate(Class = gsub("^PC P-.*", "PC-P", Class)) %>%
  mutate(Class = gsub("^PC .*", "PC", Class)) %>%
  mutate(Class = gsub("^PE O-.*", "PE-O", Class)) %>%
  mutate(Class = gsub("^PE P-.*", "PE-P", Class)) %>%
  mutate(Class = gsub("^PE .*", "PE", Class)) %>%
  mutate(Class = gsub("^PG .*", "PG", Class)) %>%
  mutate(Class = gsub("^PI .*", "PI", Class)) %>%
  mutate(Class = gsub("^PS .*", "PS", Class)) %>%
  mutate(Class = gsub("^SM.*", "SM", Class)) %>%
  mutate(Class = gsub("^DAG.*", "DAG", Class)) %>%
  mutate(Class = gsub("^TAG.*", "TAG", Class)) %>%
  mutate(Class = factor(
    Class,
    levels = c(
      "TAG",
      "SM",
      "PS",
      "PI",
      "PG",
      "PE-P",
      "PE-O",
      "PE",
      "PC-P",
      "PC-O",
      "PC",
      "LPI",
      "LPE",
      "LPC-O",
      "LPC",
      "Hex3Cer",
      "Hex2Cer",
      "HexCer",
      "DAG",
      "Cer",
      "ChE"
    )
  )) %>%
  mutate(Consensus = Consensus / 1e6, Uncertainty = Uncertainty / 1e6) # change to mol from micromol

# summarise concentration per lipid species
concentrationByLipidSpecies <- concentrationData %>%
  group_by(Lipid, Class) %>%
  summarise(MeanConsensus = mean(Consensus, na.rm = TRUE))

# sum mean consensus for each lipid class
concentrationByLipidClass <- concentrationByLipidSpecies %>%
  group_by(Class) %>%
  summarise(SumConsensus = sum(MeanConsensus, na.rm = TRUE))

#write_csv(x = concentrationByLipidSpecies, path = "plasma-mean-concentrations-lipid-species.csv")
#write_csv(x = concentrationByLipidClass, path = "plasma-summed-mean-concentration-lipid-classes.csv")

# create a combined data frame for easier plotting, distinguished by label
concentrationRangeDataByLipidSpecies <-
  concentrationByLipidSpecies %>%
  ungroup() %>%
  select(-Lipid) %>%
  rename(Concentration = MeanConsensus) %>%
  mutate(label = "Mean amount per lipid species", xend = Concentration + (0.075 * Concentration))

concentrationRangeData <- bind_rows(
  concentrationRangeDataByLipidSpecies,
  concentrationByLipidClass %>%
    ungroup() %>%
    rename(Concentration = SumConsensus) %>%
    mutate(label = "Total amount per lipid class", xend = Concentration + (0.1 * Concentration))
) %>%
  mutate(label = factor(
    label,
    levels = c("Mean amount per lipid species", "Total amount per lipid class")
  ))

#write_csv(x = concentrationRangeData, path = "plasma-concentration-range-data.csv")

# plot the mean concentrations of each lipid species and the summed concentration per lipid class
concentrationRangePlot <-
  ggplot(aes(x = Concentration), data = concentrationRangeData) +
  geom_segment(
    aes(
      x = Concentration,
      xend = xend,
      y = Class,
      yend = Class,
      color = label
    ),
    size = 4,
    show.legend = TRUE
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x)
      10 ^ x, n = 10),
    minor_breaks = log10_minor_break(),
    labels = scales::trans_format("log10", scales::math_format(10 ^ .x)),
    limits = c(1e-09, 1e-02)
  ) + annotation_logticks(sides = "bt") +
  scale_y_discrete(position = "right") +
  scale_color_manual(
    name = "",
    values = c(
      "Total amount per lipid class" = "#d95f02",
      "Mean amount per lipid species" = "#7570b3"
    )
  ) +
  ylab("") +
  xlab(bquote("Concentration" ~ group("(",paste(mol, ~L^-1), ")"))) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_rect(
      fill = "transparent",
      size = 0.5,
      linetype = "solid",
      colour = "transparent"
    ),
    legend.position = c(0.85, 0.3),
    text = element_text(size=12),
    axis.text.y = element_text(face="bold"),
    axis.title.x = element_text(face="bold"),
  )
ggsave(
  "plasmaConcentrationRangePlot.pdf",
  concentrationRangePlot,
  width = 11.69,
  height = 5.46,#8.27,
  units = "in"
)
ggsave(
  "plasmaConcentrationRangePlot.svg",
  concentrationRangePlot,
  width = 11.69,
  height = 5.46,#8.27,
  units = "in"
)

#
# scatter plots with linear fit
#

# read complete Bowden MEDM values
bowdenMedmComplete <- read_csv("Bowden-MEDM-all-species.csv") %>%
  mutate(Lipid = gsub("CE ", "ChE ", Lipid)) %>%
  mutate(Lipid = gsub("d18:0", "18:0;2", Lipid)) %>%
  mutate(Lipid = gsub("d18:1", "18:1;2", Lipid))

#unqLipids <- unique(bind_rows(qtys_nist_renamed, qtys_summary_renamed)$Lipid)
#write_tsv(data.frame(lipids=unqLipids), path="lc-lipids.tsv")
#write_tsv(data.frame(lipids=bowdenMedmComplete$Lipid), path="bowden-lipids.tsv")

# map lipidCreator names to Bowden study names
qtys_nist_renamed_trans <- qtys_nist_renamed %>%
  mutate(Lipid = gsub("Cer 18:0;2/22:0","Cer 40:0", Lipid)) %>%
  mutate(Lipid = gsub("Cer 18:0;2/24:1","Cer 42:1", Lipid)) %>%
  mutate(Lipid = gsub("Cer 18:1;2/16:0","Cer 34:1", Lipid)) %>%
  mutate(Lipid = gsub("Cer 18:1;2/18:0","Cer 36:1", Lipid)) %>%
  mutate(Lipid = gsub("Cer 18:1;2/20:0","Cer 38:1", Lipid)) %>%
  mutate(Lipid = gsub("Cer 18:1;2/22:0","Cer 40:1", Lipid)) %>%
  mutate(Lipid = gsub("Cer 18:1;2/24:0","Cer 42:1", Lipid)) %>%
  mutate(Lipid = gsub("Cer 18:1;2/24:1","Cer 42:2", Lipid)) %>%
  mutate(Lipid = gsub("DAG 14:0_18:2","DAG 32:2", Lipid)) %>%
  mutate(Lipid = gsub("DAG 16:0_18:1","DAG 34:1", Lipid)) %>%
  mutate(Lipid = gsub("DAG 16:0_18:2","DAG 34:2", Lipid)) %>% 
  mutate(Lipid = gsub("DAG 16:0_20:4","DAG 36:4", Lipid)) %>%
  mutate(Lipid = gsub("DAG 16:0_22:6","DAG 38:6", Lipid)) %>%
  mutate(Lipid = gsub("DAG 16:1_18:1","DAG 34:2", Lipid)) %>%
  mutate(Lipid = gsub("DAG 18:0_18:2","DAG 36:2", Lipid)) %>%
  mutate(Lipid = gsub("DAG 18:0_20:4","DAG 38:4", Lipid)) %>%
  mutate(Lipid = gsub("DAG 18:1_18:2","DAG 36:3", Lipid)) %>%
  mutate(Lipid = gsub("DAG 18:1_18:3","DAG 36:4", Lipid)) %>%
  mutate(Lipid = gsub("DAG 18:1_20:3","DAG 38:4", Lipid)) %>%
  mutate(Lipid = gsub("DAG 18:2_18:2","DAG 36:4", Lipid)) %>%
  mutate(Lipid = gsub("DAG 18:2_20:4","DAG 38:6", Lipid)) %>%
  mutate(Lipid = gsub("Hex2Cer 18:1;2/16:0","Hex2Cer 34:1", Lipid)) %>%
  mutate(Lipid = gsub("Hex2Cer 18:1;2/24:1","Hex2Cer 42:2", Lipid)) %>%
  mutate(Lipid = gsub("HexCer 18:1;2/16:0","HexCer 34:1", Lipid)) %>%
  mutate(Lipid = gsub("HexCer 18:1;2/18:0","HexCer 36:1", Lipid)) %>%
  mutate(Lipid = gsub("HexCer 18:1;2/20:0","HexCer 38:1", Lipid)) %>%
  mutate(Lipid = gsub("HexCer 18:1;2/22:0","HexCer 40:1", Lipid)) %>%
  mutate(Lipid = gsub("HexCer 18:1;2/24:1","HexCer 42:2", Lipid)) %>%
  mutate(Lipid = gsub("PE P-16:0/18:1","PE P-34:1", Lipid)) %>%
  mutate(Lipid = gsub("PE P-16:0/18:2","PE P-34:2", Lipid)) %>%
  mutate(Lipid = gsub("PE P-16:0/20:4","PE P-36:4", Lipid)) %>%
  mutate(Lipid = gsub("PE P-16:0/22:5","PE P-38:5", Lipid)) %>%
  mutate(Lipid = gsub("PE P-16:0/22:6","PE P-38:6", Lipid)) %>%
  mutate(Lipid = gsub("PE P-18:0/18:1","PE P-36:1", Lipid)) %>%
  mutate(Lipid = gsub("PE P-18:0/18:2","PE P-36:2", Lipid)) %>%
  mutate(Lipid = gsub("PE P-18:0/20:4","PE P-38:4", Lipid)) %>%
  mutate(Lipid = gsub("PE P-18:0/22:5","PE P-40:5", Lipid)) %>%
  mutate(Lipid = gsub("PE P-18:0/22:6","PE P-40:6", Lipid)) %>%
  mutate(Lipid = gsub("PE P-20:0/20:4","PE P-40:4", Lipid))

qtys_summary_renamed_trans <- qtys_summary_renamed %>%
  mutate(Lipid = gsub("Cer 18:0;2/22:0","Cer 40:0", Lipid)) %>%
  mutate(Lipid = gsub("Cer 18:0;2/24:1","Cer 42:1", Lipid)) %>%
  mutate(Lipid = gsub("Cer 18:1;2/16:0","Cer 34:1", Lipid)) %>%
  mutate(Lipid = gsub("Cer 18:1;2/18:0","Cer 36:1", Lipid)) %>%
  mutate(Lipid = gsub("Cer 18:1;2/20:0","Cer 38:1", Lipid)) %>%
  mutate(Lipid = gsub("Cer 18:1;2/22:0","Cer 40:1", Lipid)) %>%
  mutate(Lipid = gsub("Cer 18:1;2/24:0","Cer 42:1", Lipid)) %>%
  mutate(Lipid = gsub("Cer 18:1;2/24:1","Cer 42:2", Lipid)) %>%
  mutate(Lipid = gsub("DAG 14:0_18:2","DAG 32:2", Lipid)) %>%
  mutate(Lipid = gsub("DAG 16:0_18:1","DAG 34:1", Lipid)) %>%
  mutate(Lipid = gsub("DAG 16:0_18:2","DAG 34:2", Lipid)) %>% 
  mutate(Lipid = gsub("DAG 16:0_20:4","DAG 36:4", Lipid)) %>%
  mutate(Lipid = gsub("DAG 16:0_22:6","DAG 38:6", Lipid)) %>%
  mutate(Lipid = gsub("DAG 16:1_18:1","DAG 34:2", Lipid)) %>%
  mutate(Lipid = gsub("DAG 18:0_18:2","DAG 36:2", Lipid)) %>%
  mutate(Lipid = gsub("DAG 18:0_20:4","DAG 38:4", Lipid)) %>%
  mutate(Lipid = gsub("DAG 18:1_18:2","DAG 36:3", Lipid)) %>%
  mutate(Lipid = gsub("DAG 18:1_18:3","DAG 36:4", Lipid)) %>%
  mutate(Lipid = gsub("DAG 18:1_20:3","DAG 38:4", Lipid)) %>%
  mutate(Lipid = gsub("DAG 18:2_18:2","DAG 36:4", Lipid)) %>%
  mutate(Lipid = gsub("DAG 18:2_20:4","DAG 38:6", Lipid)) %>%
  mutate(Lipid = gsub("Hex2Cer 18:1;2/16:0","Hex2Cer 34:1", Lipid)) %>%
  mutate(Lipid = gsub("Hex2Cer 18:1;2/24:1","Hex2Cer 42:2", Lipid)) %>%
  mutate(Lipid = gsub("HexCer 18:1;2/16:0","HexCer 34:1", Lipid)) %>%
  mutate(Lipid = gsub("HexCer 18:1;2/18:0","HexCer 36:1", Lipid)) %>%
  mutate(Lipid = gsub("HexCer 18:1;2/20:0","HexCer 38:1", Lipid)) %>%
  mutate(Lipid = gsub("HexCer 18:1;2/22:0","HexCer 40:1", Lipid)) %>%
  mutate(Lipid = gsub("HexCer 18:1;2/24:1","HexCer 42:2", Lipid)) %>%
  mutate(Lipid = gsub("PE P-16:0/18:1","PE P-34:1", Lipid)) %>%
  mutate(Lipid = gsub("PE P-16:0/18:2","PE P-34:2", Lipid)) %>%
  mutate(Lipid = gsub("PE P-16:0/20:4","PE P-36:4", Lipid)) %>%
  mutate(Lipid = gsub("PE P-16:0/22:5","PE P-38:5", Lipid)) %>%
  mutate(Lipid = gsub("PE P-16:0/22:6","PE P-38:6", Lipid)) %>%
  mutate(Lipid = gsub("PE P-18:0/18:1","PE P-36:1", Lipid)) %>%
  mutate(Lipid = gsub("PE P-18:0/18:2","PE P-36:2", Lipid)) %>%
  mutate(Lipid = gsub("PE P-18:0/20:4","PE P-38:4", Lipid)) %>%
  mutate(Lipid = gsub("PE P-18:0/22:5","PE P-40:5", Lipid)) %>%
  mutate(Lipid = gsub("PE P-18:0/22:6","PE P-40:6", Lipid)) %>%
  mutate(Lipid = gsub("PE P-20:0/20:4","PE P-40:4", Lipid))

# create joined tables
nistVsMedmData <-
  left_join(bowdenMedmComplete, qtys_nist_renamed_trans, by = "Lipid")
nrow(nistVsMedmData[!is.na(nistVsMedmData$`MEDM-consensus`) & !is.na(nistVsMedmData$Consensus),])
plasmaVsMedmData <-
  left_join(bowdenMedmComplete, qtys_summary_renamed_trans, by = "Lipid")
nrow(plasmaVsMedmData[!is.na(plasmaVsMedmData$`MEDM-consensus`) & !is.na(plasmaVsMedmData$Consensus),])
nistVsSubjectsData <-
  left_join(qtys_nist_renamed, qtys_summary_renamed_trans, by = "Lipid")
nrow(nistVsSubjectsData[!is.na(nistVsSubjectsData$`Consensus.x`) & !is.na(nistVsSubjectsData$`Consensus.y`),])
# NIST SRM 1950 mean vs Bowden MEDM
# nistVsMedmData
linFormula <- y ~ x
# linear correlation between LipidCreator measurements and Bowden MEDM ring trial results
nistVsMedmPlot <-
  ggplot(nistVsMedmData, aes(x = Consensus, y = `MEDM-consensus`)) +
  geom_smooth(
    method = "lm",
    data = nistVsMedmData,
    aes(x = Consensus, y = `MEDM-consensus`),
    color = "blue",
    formula = linFormula
  ) +
  stat_poly_eq(formula = linFormula,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  geom_point(size = 3, color = "black", shape = 16) +
  theme_bw() +
  xlab(bquote("NIST SRM 1950 Mean Concentration" ~ group("(",paste(mu, mol, ~L^-1), ")"))) +
  ylab(bquote("Bowden MEDM Consensus Concentration" ~ group("(",paste(mu, mol, ~L^-1), ")"))) +
  expand_limits(y = 0) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x)
      10 ^ x, n = 4),
    minor_breaks = log10_minor_break(),
    labels = scales::trans_format("log10", scales::math_format(10 ^ .x)),
    limits = c(0.001, 10e3)
  ) + annotation_logticks(sides = "bl") +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x)
      10 ^ x, n = 4),
    minor_breaks = log10_minor_break(),
    labels = scales::trans_format("log10", scales::math_format(10 ^ .x)),
    limits = c(0.001, 10e3)
  )
ggsave(
  "nistVsMedmPlotLinearCorrelation.pdf",
  nistVsMedmPlot,
  width = 11.69,
  height = 8.27,
  units = "in"
)
ggsave(
  "nistVsMedmPlotLinearCorrelation.svg",
  nistVsMedmPlot,
  width = 11.69,
  height = 8.27,
  units = "in"
)
ggsave(
  "nistVsMedmPlotLinearCorrelation.png",
  nistVsMedmPlot,
  width = 11.69,
  height = 8.27,
  units = "in"
)

# NIST SRM 1950 mean vs LipidCreator Mean
# nistVsSubjectsData
nistVsSubjectsPlot <-
  ggplot(nistVsSubjectsData, aes(x = `Consensus.x`, y = `Consensus.y`)) +
  geom_smooth(
    method = "lm",
    data = nistVsSubjectsData,
    aes(x = `Consensus.x`, y = `Consensus.y`),
    color = "blue",
    formula = linFormula
  ) +
  stat_poly_eq(formula = linFormula,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  geom_point(size = 3, color = "black", shape = 16) +
  theme_bw() +
  xlab(bquote("NIST SRM 1950 Mean Concentration" ~ group("(",paste(mu, mol, ~L^-1), ")"))) +
  ylab(bquote("Human Plasma Mean Concentration" ~ group("(",paste(mu, mol, ~L^-1), ")"))) +
  expand_limits(y = 0) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x)
      10 ^ x, n = 4),
    minor_breaks = log10_minor_break(),
    labels = scales::trans_format("log10", scales::math_format(10 ^ .x)),
    limits = c(0.001, 10e3)
  ) + annotation_logticks(sides = "bl") +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x)
      10 ^ x, n = 4),
    minor_breaks = log10_minor_break(),
    labels = scales::trans_format("log10", scales::math_format(10 ^ .x)),
    limits = c(0.001, 10e3)
  )
ggsave(
  "nistVsSubjectsPlotLinearCorrelation.pdf",
  nistVsSubjectsPlot,
  width = 11.69,
  height = 8.27,
  units = "in"
)
ggsave(
  "nistVsSubjectsPlotLinearCorrelation.svg",
  nistVsSubjectsPlot,
  width = 11.69,
  height = 8.27,
  units = "in"
)
ggsave(
  "nistVsSubjectsPlotLinearCorrelation.png",
  nistVsSubjectsPlot,
  width = 11.69,
  height = 8.27,
  units = "in"
)

# Plasma Mean vs Bowden MEDM
# plasmaVsMedmData
plasmaVsMedmPlot <-
  ggplot(plasmaVsMedmData, aes(x = Consensus, y = `MEDM-consensus`)) +
  geom_smooth(
    method = "lm",
    data = plasmaVsMedmData,
    aes(x = Consensus, y = `MEDM-consensus`),
    color = "blue",
    formula = linFormula
  ) +
  stat_poly_eq(formula = linFormula,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  geom_point(size = 3, color = "black", shape = 16) +
  theme_bw() +
  xlab(bquote("Human Plasma Mean Concentration" ~ group("(",paste(mu, mol, ~L^-1), ")"))) +
  ylab(bquote("Bowden MEDM Consensus Concentration" ~ group("(",paste(mu, mol, ~L^-1), ")"))) +
  expand_limits(y = 0) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x)
      10 ^ x, n = 4),
    minor_breaks = log10_minor_break(),
    labels = scales::trans_format("log10", scales::math_format(10 ^ .x)),
    limits = c(0.001, 10e3)
  ) + annotation_logticks(sides = "bl") +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x)
      10 ^ x, n = 4),
    minor_breaks = log10_minor_break(),
    labels = scales::trans_format("log10", scales::math_format(10 ^ .x)),
    limits = c(0.001, 10e3)
  )
ggsave(
  "plasmaVsMedmPlotLinearCorrelation.pdf",
  plasmaVsMedmPlot,
  width = 11.69,
  height = 8.27,
  units = "in"
)
ggsave(
  "plasmaVsMedmPlotLinearCorrelation.svg",
  plasmaVsMedmPlot,
  width = 11.69,
  height = 8.27,
  units = "in"
)
ggsave(
  "plasmaVsMedmPlotLinearCorrelation.png",
  plasmaVsMedmPlot,
  width = 11.69,
  height = 8.27,
  units = "in"
)