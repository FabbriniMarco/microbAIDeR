data <- read.delim("inst/example_data/genera_example_2.csv", header=T, row.names=1, sep=',')
env <- read.delim("inst/example_data/env_example.csv", header=T, row.names = 1, sep=',')

library(microbAIDeR)

group = factor(env$Group, levels=c("NEG", "POS"))
group = c("Group", "Timepoint")

env[,"Group"] = factor(env[,"Group"], levels=c("NEG", "POS"))
env[,"Timepoint"] = factor(env[,"Timepoint"], levels=c("T0", "T1"))

# compute_LMM_and_plot(
#   data = data,
#   group = group,
#   taxlevel = "test",
#   save.path = "Test_outs",
#   color.grouping = c("firebrick3", "goldenrod2", "dodgerblue3", "forestgreen"),
#   mode = "taxa",
#   meta = env,
#   fixed_effects = c("Param.1", "Param.2", "Group"),
#   random_effects = c("Param.3", "Param.4")
# )



taxlevel = "test"
save.path = "Test_outs"
color.grouping = c("firebrick3", "goldenrod2", "dodgerblue3", "forestgreen")
comparison.list = list( c("NEG T0", "NEG T1"), c("NEG T0", "POS T0"), c("POS T0", "POS T1"), c("NEG T1", "POS T1") )
p.adjust.method = "fdr"
trends = TRUE
plot.not.sig = TRUE
paired = FALSE
mode = c("wholetable")     # singlefeature vs wholetable
data_type = "relabb"
meta = env                  # data.frame with covariates; rownames = sample ids (must match colnames(data))
fixed_effects = c("Group * Timepoint * taxa_name", "age", "gender")
random_effects = NULL
dispformula = TRUE                   # TRUE => ~ 0 + taxa_var (if provided) else ~ 1; FALSE/NULL => omit; character => custom
ziformula = TRUE                    # TRUE => ~ 0 + taxa_var (if provided) else ~ 1; FALSE/NULL => omit; character => custom
family_model = glmmTMB::beta_family(link = "logit")
var.timepoint.random.slope = "Timepoint"    # e.g., "New_timepoint"
var.subject.random.intercept = "Subject_ID" # e.g., "Subject_ID"
include_random_intercept_by_taxa = TRUE
conflevel = 0.95
save_single_model_summaries = TRUE
nrow.graph = 2
ncol.graph = 2
width.graph = 4.5
height.graph = 3.5
horiz = FALSE
ggplot.margins = c(.18, .18, .18, .6)
box.lwd = 0.4
jitter.pch = 21
jitter.stroke = 0.15
jitter.size = 0.7
jitter.color = "grey22"
signif.step.increase = 0.12
signif.text.size = 3
signif.line.size = 0.4
contrast.color = "ivory1"
text.x.size = 6
text.y.size = 6
text.y.title.size = 8
smoothing = FALSE
smoothing.lwd = 1
smoothing.color = "darkred"
smoothing.se = FALSE
smoothing.method="loess"
additional.params = NULL
align.legend = FALSE
plot.order = "kruskal"
pattern.fill = FALSE
pattern = "stripe"
pattern.angle = 45
pattern.alpha = 0.4
pattern.density = 0.1
pattern.spacing = 0.05





