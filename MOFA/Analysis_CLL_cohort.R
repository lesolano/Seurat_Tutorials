library(MOFA2)
library(data.table)
library(randomForest)
library(psych)
library(tidyverse)
# (Optional) set up reticulate connection with Python
# library(reticulate)
# reticulate::use_python("/Users/ricard/anaconda3/envs/base_new/bin/python", required = T)

###############
## Load data ##
###############

# Multiple formats are allowed for the input data:

## -- Option 1 -- ##
# nested list of matrices, where the first index refers to the view and the second index refers to the group.
# samples are stored in the rows and features are stored in the columns.
# Missing values must be filled with NAs, including samples missing an entire view

# (...)

## -- Option 2 -- ##
# data.frame with columns ["sample","feature","view","group","value"]
# In this case there is no need to have missing values in the data.frame,
# they will be automatically filled in when creating the corresponding matrices

file = "ftp://ftp.ebi.ac.uk/pub/databases/mofa/getting_started/data.txt.gz"
data = fread(file)

#######################
# Create MOFA object ##
#######################

MOFAobject <- create_mofa(data)

# Visualise data structure
plot_data_overview(MOFAobject)

####################
## Define options ##
####################

# Data options
# - scale_views: if views have very different ranges/variances, it is good practice to scale each view to unit variance (default is FALSE)
data_opts <- get_default_data_options(MOFAobject)


# Model options
# - likelihoods: likelihood per view (options are "gaussian","poisson","bernoulli"). "gaussian" is used by default
# - num_factors: number of factors. By default K=10
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10

# Training options
# - maxiter: number of iterations
# - convergence_mode: "fast", "medium", "slow". For exploration, the fast mode is good enough.
# - drop_factor_threshold: minimum variance explained criteria to drop factors while training. Default is -1 (no dropping of factors)
# - gpu_mode: use GPU mode? This needs cupy installed and a functional GPU, see https://biofam.github.io/MOFA2/gpu_training.html
# - seed: random seed
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "medium"
train_opts$seed <- 42

#########################
## Prepare MOFA object ##
#########################

MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)

#####################
## Train the model ##
#####################

MOFAobject <- run_mofa(MOFAobject, outfile = paste0(getwd(),"/model.hdf5"))

####################
## Not confident the trained model is generating expected output ##
####################

MOFAobject <- H5Fopen('model.hdf5')

# Load precomputed model
MOFAobject <- readRDS(url("http://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/MOFA2_CLL.rds"))

#Vignette steps
slotNames(MOFAobject)
names(MOFAobject@data)
dim(MOFAobject@data$Drugs$group1)
names(MOFAobject@expectations)
# Dimensionality of the factor matrix: 200 samples, 15 factors
dim(MOFAobject@expectations$Z$group1)
# Dimensionality of the mRNA Weight matrix: 5000 features, 15 factors
dim(MOFAobject@expectations$W$mRNA)

CLL_metadata <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt")
# Sanity check
# stopifnot(all(sort(CLL_metadata$sample)==sort(unlist(samples_names(MOFAobject)))))

# Add sample metadata to the model
samples_metadata(MOFAobject) <- CLL_metadata

plot_factor_cor(MOFAobject)

plot_variance_explained(MOFAobject, max_r2=15)

plot_variance_explained(MOFAobject, plot_total = T)[[2]]

correlate_factors_with_covariates(MOFAobject, 
                                  covariates = c("Gender","died","age"), 
                                  plot="log_pval"
)

#Characterization of Factor 1
correlate_factors_with_covariates(MOFAobject, 
                                  covariates = c("Gender","died","age"), 
                                  plot="log_pval"
)

plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Factor1"
)

plot_weights(MOFAobject,
             view = "Mutations",
             factor = 1,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

plot_top_weights(MOFAobject,
                 view = "Mutations",
                 factor = 1,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)

plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "IGHV",
            add_violin = TRUE,
            dodge = TRUE
)

plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Gender",
            dodge = TRUE,
            add_violin = TRUE
)

plot_weights(MOFAobject, 
             view = "mRNA", 
             factor = 1, 
             nfeatures = 10
)

plot_data_scatter(MOFAobject, 
                  view = "mRNA",
                  factor = 1,  
                  features = 4,
                  sign = "positive",
                  color_by = "IGHV"
) + ggplot2::labs(y="RNA expression")

plot_data_scatter(MOFAobject, 
                  view = "mRNA",
                  factor = 1,  
                  features = 4,
                  sign = "negative",
                  color_by = "IGHV"
) + ggplot2::labs(y="RNA expression")

plot_data_heatmap(MOFAobject, 
                  view = "mRNA",
                  factor = 1,  
                  features = 25,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)

plot_data_heatmap(MOFAobject, 
                  view = "mRNA",
                  factor = 1,  
                  features = 25,
                  denoise = TRUE,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)

#Characterization of Factor 3
plot_weights(MOFAobject, 
             view = "Mutations", 
             factor = 3, 
             nfeatures = 10,
             abs = F
)

plot_factor(MOFAobject, 
            factors = 3, 
            color_by = "trisomy12",
            dodge = TRUE,
            add_violin = TRUE
)

plot_data_scatter(MOFAobject, 
                  view = "Drugs",
                  factor = 3,  
                  features = 4,
                  sign = "positive",
                  color_by = "trisomy12"
) + ggplot2::labs(y="Drug response (cell viability)")

plot_data_heatmap(MOFAobject, 
                  view = "mRNA",
                  factor = 3,  
                  features = 25,
                  denoise = TRUE,
                  cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)

#Inspection of combinations of Factors
p <- plot_factors(MOFAobject, 
                  factors = c(1,3), 
                  color_by = "IGHV",
                  shape_by = "trisomy12",
                  dot_size = 2.5,
                  show_missing = T
)

p <- p + 
  ggplot2::geom_hline(yintercept=-1, linetype="dashed") +
  ggplot2::geom_vline(xintercept=(-0.5), linetype="dashed")

print(p)

#Prediction of clinical subgroups
suppressPackageStartupMessages(library(randomForest))

# Prepare data
df <- as.data.frame(get_factors(MOFAobject, factors=c(1,2))[[1]])

# Train the model for IGHV
df$IGHV <- as.factor(MOFAobject@samples_metadata$IGHV)
model.ighv <- randomForest(IGHV ~ ., data=df[!is.na(df$IGHV),], ntree=10)
df$IGHV <- NULL

# Do predictions
MOFAobject@samples_metadata$IGHV.pred <- stats::predict(model.ighv, df)

# Train the model for Trisomy12
df$trisomy12 <- as.factor(MOFAobject@samples_metadata$trisomy12)
model.trisomy12 <- randomForest(trisomy12 ~ ., data=df[!is.na(df$trisomy12),], ntree=10)
df$trisomy12 <- NULL

MOFAobject@samples_metadata$trisomy12.pred <- stats::predict(model.trisomy12, df)

#Plot predictions for IGHV
MOFAobject@samples_metadata$IGHV.pred_logical <- c("True","Predicted")[as.numeric(is.na(MOFAobject@samples_metadata$IGHV))+1]

p <- plot_factors(MOFAobject, 
                  factors = c(1,3), 
                  color_by = "IGHV.pred",
                  shape_by = "IGHV.pred_logical",
                  dot_size = 2.5,
                  show_missing = T
)

p <- p + 
  ggplot2::geom_hline(yintercept=-1, linetype="dashed") +
  ggplot2::geom_vline(xintercept=(-0.5), linetype="dashed")

print(p)

