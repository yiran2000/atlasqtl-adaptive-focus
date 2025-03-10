library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(ggplot2)
library(atlasqtl)
library(tictoc)
library(echoseq)
library(PRROC)
library(patchwork)
library(stringr)


# simulation_tool_path = "~/project/response_simulation"
simulation_tool_path = "C:/Users/Yiran/OneDrive - University of Cambridge/Documents/PhD/response_simulation"
source(file.path(simulation_tool_path, "data_simulation_tools/simulate_withRealGenome.R"))
source("reshape_atlasQTL_res.R")

################################################################################
#read in pre-save data
#X
list = readRDS(file.path(simulation_tool_path, "data/sim_list_3.rds"))
X_real = as.matrix(list$X)

#Y_real_corr
Y_real_corr = readRDS(file.path(simulation_tool_path, "data/Y_real_corr.rds"))

#protein rank
# protein_rank = readLines("data/protein_rank.txt")

#snp gene info
snp_gene_info = fread(file.path(simulation_tool_path, "data/snp_gene_info.csv"))

################################################################################
#simulate one single hotspot
list_sim = simulate_withRealGenome(
  X = X_real[,801:1000], 
  q = 1000,
  protein_ls = paste0("protein_", 1:1000),
  active_protein = NULL,
  active_ratio_p = 2/200, 
  active_ratio_q = 0.05,
  residual_cor_mat = NULL,
  vec_rho_phenos = NULL,
  missing_ratio = 0, 
  max_tot_pve = NULL, 
  nb_phenos_per_block =  200, 
  min_cor = 0, max_cor = 0.9,
  m = 0.02,
  sh2 = 5,
  seed  = 2024
)


X_sim = list_sim$snps
Y_sim = list_sim$phenos
beta_sim = list_sim$beta
pat_sim = list_sim$pat

# Y_sim_corr = cor(Y_sim)
# hist(as.vector(beta_sim)[as.vector(beta_sim)!=0])

#visualize simulated associations
data.frame(ID = rownames(pat_sim), ProteinCount = rowSums(pat_sim)) %>% 
  left_join(snp_gene_info, by = "ID") %>% 
  ggplot(
    aes(x = POS, y = ProteinCount)
  ) +
  geom_point()+
  theme_bw()+
  ggtitle("Simulated # of associated proteins by SNP")


################################################################################
#run atlasQTL

# p = ncol(X_sim)
# q = 2919
# n = nrow(X_sim)


mu_t = 1
v_t = 4

devtools::load_all()
system.time(res_atlas =  atlasqtl(Y = as.matrix(Y_sim), X = as.matrix(X_sim),
                         p0 = c(mu_t, v_t),
                         user_seed = 1, maxit=1000,
                         batch = "y",
                         tol_init = 0.1,
                         tol_loose = 0.1,
                         tol_tight = 0.1,
                         burn_in = 10,
                         maxit_full = 5,
                         maxit_subsample = 100000,
                         n_partial_update = 500,
                         eval_perform = T,
                         epsilon= c(2, 1.5, 0.25),
                         partial_elbo = F, #whether we calculate elbo with all responses or only a part
                         partial_elbo_eval = F, #whether diff_lb = lb_new -lb_old or (lb_new-lb_old)/length(sample_q)
                         thinned_elbo_eval = T))

#evaluate AUROC and AUPRC
roc <- PRROC::roc.curve(scores.class0 = as.vector(res_atlas$gam_vb), weights.class0 = as.vector(as.matrix(pat_sim)), curve =TRUE)
AUROC = roc$auc
AUROC 

pr <- PRROC::pr.curve(scores.class0 = as.vector(res_atlas$gam_vb), weights.class0 = as.vector(as.matrix(pat_sim)), curve =TRUE)
AUPRC = pr$auc.integral
AUPRC


################################################################################

# manhattan plot
# simulated true associations
p1 = data.frame(ID = rownames(pat_sim), ProteinCount = rowSums(pat_sim)) %>% 
  left_join(snp_gene_info, by = "ID") %>% 
  ggplot(
    aes(x = POS, y = ProteinCount)
  ) +
  # scale_x_continuous(limits = gene_range)+
  geom_point()+
  theme_bw()

# atlasQTL inferred associations
res_simData = reshape_atlasQTL_res(obj_atlasqtl,snp_gene_info = snp_gene_info)
p2 = res_simData[[1]]%>% 
  group_by(POS, ID) %>% 
  summarise(n_assoc = sum(corr_metric > 0.9)) %>% 
  ggplot(aes(x = POS, y = n_assoc))+
  geom_point()+
  theme_bw()

p3 = res_simData[[2]] %>% 
  ggplot(aes(x = POS, y = theta))+
  geom_point()+
  theme_bw()

wrap_plots(p1, p2, p3, ncol = 1)

################################################################################
#FDR
# exact FDR
gam_vb = obj_atlasqtl$gam_vb
predicted <- gam_vb > 0.9

TP <- sum((predicted==1) & (pat_sim == 1))
FP <- sum(predicted & (pat_sim == 0))
FDR <- FP / sum(predicted==1)
FDR

power = TP/sum(pat_sim == 1)
power


# loci-wise FDR
active_set_pred = which(apply(gam_vb, 2, function(col) any(col > 0.9)))
active_set_true = which(apply(pat_sim, 2, function(col) any(col == 1)))
TP_proteins = sum(active_set_pred %in% active_set_true)
FP_proteins = sum(!active_set_pred %in% active_set_true)

FP_proteins / (FP_proteins + TP_proteins)

TP_proteins/1000

