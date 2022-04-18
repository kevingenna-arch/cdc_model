### File to generate plots for the documents behind june's conf


# 0 - settings and run the matlab code -------------------------------------------------

library(tidyverse)
library(R.matlab)
library(matlabr)
library(xtable)
options(xtable.include.rownames = F,
        xtable.floating = F,
        xtable.latex.environment = F)

# find matlab and run the script generating the files
get_matlab()
run_matlab_script("../cdc_juin/tables_rondes.m",
                  verbose = F)
yrs <- c(2020, 2030, 2050, 2070, 2100)


# TR0: some background ----------------------------------------------------

# effect of productivity

tech_dyn <- read_csv("../output/tr0_tech.csv",
                     show_col_types = FALSE) %>%
  separate(col = variable,
           into = c("var", "skill", "age"),
           sep = '_',
           remove = F) %>%
  mutate(age = as.numeric(age)*5+15,
         period = period*5 + 1900,
         mode = case_when(
           mod == 'base' ~ "dA=1.3%",
           mod == 'opt' ~ "dA=1.8%",
           mod == 'pess' ~ 'dA=1%',
           mod == 'hist' ~ 'Donn. Hist.',
           mod == 'low' ~ 'dA=.25%'
         )
         )

plot_tech_dyn_y <- tech_dyn %>%
  filter(var == 'y') %>%
  ggplot(aes(x = period, y = value, colour = mode)) +
  geom_line() + xlab('') + ylab('PIB') + theme_minimal() +
  theme(legend.title = element_blank(), legend.position="bottom")

ggsave(filename = "tr0_tech_dyn_pib.eps",
       plot = plot_tech_dyn_y,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)


# TR1: labour supply -----------------------------------------------------

##### SS comparison 1 ####
pens_ss <- read_csv('../output/tr1_pensions_ss.csv',
                    show_col_types = F) %>%
  pivot_wider(names_from = 'mod',
              values_from = 'value') %>%
  separate(variable,
           into = c('vars', 'skill', 'cohort'),
           sep = '_',
           remove = F,
           convert = T) %>%
  mutate(skill = factor(skill, c('NQ', 'Q')),
         cohort = cohort*5+15,
         across(.cols = c(age_55, age_65, age_70),
                .fns = ~ (.x - age_60)*100/age_60,
                .names = "{.col}_pct"))

# SS tables -- levels

pens_ss_lvl_macro <- pens_ss %>%
  filter(variable %in% c('y', 'nbar', 'c', 'L_Q', 'L_NQ', 'kbar', 'H', 'I', 'R', 'rd')) %>%
  select(!contains('_pct'), -cohort, -variable)

pens_ss_lvl_macro %>%
  xtable() %>%
  print(type = 'latex',
        file = '../output/tr1_ss_macro_lvl.tex',
        floating = F,
        latex.environment = F)

pens_ss_lvl_fiscal <- pens_ss %>%
  filter(grepl('ratio', vars) |
           grepl('tau', vars) |
           grepl('T', vars) |
           vars %in% c("g", "D", "penbase", "penind", "retire", "Deped")) %>%
  select(!contains('_pct'), -cohort, -variable)


pens_ss_lvl_fiscal %>%
  xtable() %>%
  print(type = 'latex',
        file = '../output/tr1_ss_fiscal_lvl.tex',
        floating = F,
        latex.environment = F)

# SS tables -- pct
pens_ss_pct_macro <- pens_ss %>%
  filter(variable %in% c('y', 'nbar', 'c', 'L_Q', 'L_NQ', 'kbar', 'H', 'I', 'R', 'rd')) %>%
  select(vars, skill, ends_with('_pct'), -cohort, -variable)

pens_ss_pct_macro %>%
  xtable() %>%
  print(type = 'latex',
        file = '../output/tr1_ss_macro_pct.tex',
        floating = F,
        latex.environment = F)

pens_ss_pct_fiscal <- pens_ss %>%
  filter(grepl('ratio', vars) |
           grepl('tau', vars) |
           grepl('T', vars) |
           vars %in% c("g", "D", "penbase", "penind", "retire", "Deped")) %>%
  select(vars, skill, contains('_pct'), -cohort, -variable)


pens_ss_pct_fiscal %>%
  xtable() %>%
  print(type = 'latex',
        file = '../output/tr1_ss_fiscal_pct.tex',
        floating = F,
        latex.environment = F)

# life-cycle plots -- levels
pens_ss_long_lvl <- pens_ss %>%
  select(!ends_with('_pct')) %>%
  pivot_longer(cols = contains('age_'),
               names_to = 'model',
               values_to = 'value') %>%
  mutate(model = gsub('age_', '', model))

plot_pens_ss_welf <- pens_ss_long_lvl %>%
  filter(vars == 'welf') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_pens_ss_c <- pens_ss_long_lvl %>%
  filter(vars == 'c' & variable != 'c') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_pens_ss_montret <- pens_ss_long_lvl %>%
  filter(vars == 'pen') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_pens_ss_k <- pens_ss_long_lvl %>%
  filter(vars == 'k') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_pens_ss_wage <- pens_ss_long_lvl %>%
  filter(vars == 'w') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

ggsave(filename = "tr1_pens_ss_lvl_welf.eps",
       plot = plot_pens_ss_welf,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_pens_ss_lvl_conso.eps",
       plot = plot_pens_ss_c,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_pens_ss_lvl_montret.eps",
       plot = plot_pens_ss_montret,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_pens_ss_lvl_k.eps",
       plot = plot_pens_ss_k,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_pens_ss_lvl_wage.eps",
       plot = plot_pens_ss_wage,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

# life-cycle plots -- pct
pens_ss_long_pct <- pens_ss %>%
  select(-(5:8)) %>%
  pivot_longer(cols = contains('age_'),
               names_to = 'model',
               values_to = 'value') %>%
  mutate(model = gsub('age_', '', model),
         model = gsub('_pct', '', model),
         model = factor(model))

plot_pens_ss_pct_welf <- pens_ss_long_pct %>%
  filter(vars == 'welf') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  geom_hline(yintercept = 0) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_pens_ss_pct_c <- pens_ss_long_pct %>%
  filter(vars == 'c' & variable != 'c') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  geom_hline(yintercept = 0) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_pens_ss_pct_montret <- pens_ss_long_pct %>%
  filter(vars == 'pen') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  geom_hline(yintercept = 0) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_pens_ss_pct_k <- pens_ss_long_pct %>%
  filter(vars == 'k') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  geom_hline(yintercept = 0) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_pens_ss_pct_wage <- pens_ss_long_pct %>%
  filter(vars == 'w') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  geom_hline(yintercept = 0) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

ggsave(filename = "tr1_pens_ss_pct_welf.eps",
       plot = plot_pens_ss_pct_welf,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_pens_ss_pct_conso.eps",
       plot = plot_pens_ss_pct_c,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_pens_ss_pct_montret.eps",
       plot = plot_pens_ss_pct_montret,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_pens_ss_pct_k.eps",
       plot = plot_pens_ss_pct_k,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_pens_ss_pct_wage.eps",
       plot = plot_pens_ss_pct_wage,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)


##### SS + pens ratio comparison ####
pensdep_ss <- read_csv('../output/tr1_pensions_dep_ss.csv',
                    show_col_types = F) %>%
  pivot_wider(names_from = c('mod', 'def'),
              values_from = 'value') %>%
  rename("age_60"="age_60_0.14",
         "age_55"="age_55_0.17",
         "age_65"="age_65_0.12",
         "age_70"="age_70_0.1") %>%
  separate(variable,
           into = c('vars', 'skill', 'cohort'),
           sep = '_',
           remove = F,
           convert = T) %>%
  mutate(skill = factor(skill, c('NQ', 'Q')),
         cohort = cohort*5+15,
         across(.cols = c(age_55, age_65, age_70),
                .fns = ~ (.x - age_60)*100/age_60,
                .names = "{.col}_pct"))

##### SS+dep tables -- levels ####

pens_ssdep_lvl_macro <- pensdep_ss %>%
  filter(variable %in% c('y', 'nbar', 'c', 'L_Q', 'L_NQ', 'kbar', 'H', 'I', 'R', 'rd')) %>%
  select(!contains('_pct'), -cohort, -variable)

pens_ssdep_lvl_macro %>%
  xtable() %>%
  print(type = 'latex',
        file = '../output/tr1_ssdep_macro_lvl.tex',
        floating = F,
        latex.environment = F)

pens_ssdep_lvl_fiscal <- pensdep_ss %>%
  filter(grepl('ratio', vars) |
           grepl('tau', vars) |
           grepl('T', vars) |
           vars %in% c("g", "D", "penbase", "penind", "retire", "Deped")) %>%
  select(!contains('_pct'), -cohort, -variable)


pens_ssdep_lvl_fiscal %>%
  xtable() %>%
  print(type = 'latex',
        file = '../output/tr1_ssdep_fiscal_lvl.tex',
        floating = F,
        latex.environment = F)

##### SS+dep tables -- pct ####
pens_ssdep_pct_macro <- pensdep_ss %>%
  filter(variable %in% c('y', 'nbar', 'c', 'L_Q', 'L_NQ', 'kbar', 'H', 'I', 'R', 'rd')) %>%
  select(vars, skill, ends_with('_pct'), -cohort, -variable)

pens_ssdep_pct_macro %>%
  xtable() %>%
  print(type = 'latex',
        file = '../output/tr1_ssdep_macro_pct.tex',
        floating = F,
        latex.environment = F)

pens_ssdep_pct_fiscal <- pensdep_ss %>%
  filter(grepl('ratio', vars) |
           grepl('tau', vars) |
           grepl('T', vars) |
           vars %in% c("g", "D", "penbase", "penind", "retire", "Deped")) %>%
  select(vars, skill, contains('_pct'), -cohort, -variable)


pens_ssdep_pct_fiscal %>%
  xtable() %>%
  print(type = 'latex',
        file = '../output/tr1_ssdep_fiscal_pct.tex',
        floating = F,
        latex.environment = F)


##### life-cycle plots -- levels ####
pens_ssdep_long_lvl <- pensdep_ss %>%
  select(!ends_with('_pct')) %>%
  pivot_longer(cols = contains('age_'),
               names_to = 'model',
               values_to = 'value') %>%
  mutate(model = gsub('age_', '', model))

plot_pens_ssdep_welf <- pens_ssdep_long_lvl %>%
  filter(vars == 'welf') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_pens_ssdep_c <- pens_ssdep_long_lvl %>%
  filter(vars == 'c' & variable != 'c') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_pens_ssdep_montret <- pens_ssdep_long_lvl %>%
  filter(vars == 'pen') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_pens_ssdep_k <- pens_ssdep_long_lvl %>%
  filter(vars == 'k') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_pens_ssdep_wage <- pens_ssdep_long_lvl %>%
  filter(vars == 'w') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

ggsave(filename = "tr1_pens_ssdep_lvl_welf.eps",
       plot = plot_pens_ssdep_welf,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_pens_ssdep_lvl_conso.eps",
       plot = plot_pens_ssdep_c,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_pens_ssdep_lvl_montret.eps",
       plot = plot_pens_ssdep_montret,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_pens_ssdep_lvl_k.eps",
       plot = plot_pens_ssdep_k,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_pens_ssdep_lvl_wage.eps",
       plot = plot_pens_ssdep_wage,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

##### life-cycle plots -- pct ####
pens_ssdep_long_pct <- pensdep_ss %>%
  select(-(5:8)) %>%
  pivot_longer(cols = contains('age_'),
               names_to = 'model',
               values_to = 'value') %>%
  mutate(model = gsub('age_', '', model),
         model = gsub('_pct', '', model),
         model = factor(model))

plot_pens_ssdep_pct_welf <- pens_ssdep_long_pct %>%
  filter(vars == 'welf') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  geom_hline(yintercept = 0) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_pens_ssdep_pct_c <- pens_ssdep_long_pct %>%
  filter(vars == 'c' & variable != 'c') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  geom_hline(yintercept = 0) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_pens_ssdep_pct_montret <- pens_ssdep_long_pct %>%
  filter(vars == 'pen') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  geom_hline(yintercept = 0) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_pens_ssdep_pct_k <- pens_ssdep_long_pct %>%
  filter(vars == 'k') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  geom_hline(yintercept = 0) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_pens_ssdep_pct_wage <- pens_ssdep_long_pct %>%
  filter(vars == 'w') %>%
  ggplot(aes(x = cohort, y = value, colour = model))+
  geom_line() + facet_wrap(.~skill, nrow = 2) +
  geom_hline(yintercept = 0) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

ggsave(filename = "tr1_pens_ssdep_pct_welf.eps",
       plot = plot_pens_ssdep_pct_welf,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_pens_ssdep_pct_conso.eps",
       plot = plot_pens_ssdep_pct_c,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_pens_ssdep_pct_montret.eps",
       plot = plot_pens_ssdep_pct_montret,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_pens_ssdep_pct_k.eps",
       plot = plot_pens_ssdep_pct_k,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_pens_ssdep_pct_wage.eps",
       plot = plot_pens_ssdep_pct_wage,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)


##### Dynamics for pensions ####

pens_dyn <- read_csv('../output/tr1_pensions_dyna.csv',
                     show_col_types = F) %>%
  separate(variable,
           into = c('vars', 'skill', 'cohort'),
           sep = '_',
           convert = T,
           remove = F) %>%
  mutate(period = period*5+1900,
         mod = gsub('age_', '', mod),
         mod = factor(mod),
         skill = factor(skill),
         cohort = cohort*5+15)

pens_dyn_pib <- pens_dyn %>%
  filter(vars == 'y') %>%
  select(-c(variable, skill, cohort, vars)) %>%
  # normalise to 1 in 2020, by model
  group_by(mod) %>%
  mutate(value_norm = value/value[period==2020]) %>%
  ungroup() %>%
  # pct deviations from baseline
  mutate(value_pctdev = 100*(value_norm/value_norm[mod == 60] - 1)) %>%
  group_by(mod) %>%
  mutate(mm = mean(value_pctdev[period>=2020])) %>%
  ungroup()

plot_pens_dyn_pib_norm <- pens_dyn_pib %>%
  filter(period>=2020) %>%
  ggplot(aes(x = period, y = (value_norm - 1), colour = mod)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")



plot_pens_dyn_pib_pct <- pens_dyn_pib %>%
  filter(period>=2020) %>%
  ggplot(aes(x = period, y = (value_pctdev), colour = mod)) +
  geom_line(aes(y = mm, colour = mod), linetype = 'dashed') +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")



ggsave(filename = "tr1_pens_dyn_pib_norm.eps",
       plot = plot_pens_dyn_pib_norm,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_pens_dyn_pib_norm_pct.eps",
       plot = plot_pens_dyn_pib_pct,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

## tables

pens_dyn_pib %>%
  filter(period %in% yrs) %>%
  select(mod, period, value_pctdev) %>%
  pivot_wider(names_from = period, values_from = value_pctdev) %>%
  rename(Age = mod) %>%
  xtable() %>%
  print(type = 'latex',
        file = '../output/tr1_pens_dyn_pib_pctdev.tex')


# TR1 -- shocks migration -------------------------------------------------

mig_dyn <- read_csv('../output/tr1_labour_mig.csv',
                    show_col_types = F) %>%
  separate(variable,
           into = c("vars", "skill", "cohort"),
           sep = '_',
           convert = T) %>%
  mutate(cohort = cohort*5+15,
         period = period*5+1900,
         skill = factor(skill),
         mod = case_when(mod == 'base' ~ 'Base',
                         mod == 'shock' ~ 'Red. Chom. 10%'),
         mod = factor(mod),
         value_norm = 100*(value/value[mod == 'Base'] - 1))

plot_mig_dyn_pib <- mig_dyn %>%
  filter(vars == 'y') %>%
  ggplot(aes(x = period, y = value, colour = mod)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_mig_dyn_lab <- mig_dyn %>%
  filter(vars == 'L') %>%
  ggplot(aes(x = period,
             y = value,
             colour = mod,
             group = interaction(skill, mod),
             linetype = skill)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_mig_dyn_ptot <- mig_dyn %>%
  filter(vars == 'Ptot') %>%
  ggplot(aes(x = period, y = value, colour = mod)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_mig_dyn_h <- mig_dyn %>%
  filter(vars == 'H') %>%
  ggplot(aes(x = period, y = value, colour = mod)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_mig_dyn_nbar <- mig_dyn %>%
  filter(vars == 'nbar') %>%
  ggplot(aes(x = period, y = value, colour = mod)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")


ggsave(filename = "tr1_mig_dyn_pib.eps",
       plot = plot_mig_dyn_pib,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_mig_dyn_lab.eps",
       plot = plot_mig_dyn_lab,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_mig_dyn_pop.eps",
       plot = plot_mig_dyn_ptot,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_mig_dyn_h.eps",
       plot = plot_mig_dyn_h,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr1_mig_dyn_nbar.eps",
       plot = plot_mig_dyn_nbar,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

# tables

mig_dyn %>%
  filter(vars %in% c('nbar', 'y', 'H', 'L', 'Ptot') &
           period %in% yrs &
           mod != 'Base') %>%
  select(-cohort, -value, -mod) %>%
  pivot_wider(names_from = 'period', values_from = 'value_norm') %>%
  select(-3) %>%
  arrange(desc(vars), skill) %>%
  xtable() %>%
  print(type = 'latex', file = '../output/tr1_mig_dyn.tex')



# TR2: competences --------------------------------------------------------


prodind_ss <- read_csv('../output/tr2_skill_ss.csv',
                       show_col_types = F) %>%
  separate(variable,
           into = c("vars", "skill", "cohort"),
           sep = '_',
           remove = F,
           convert = T) %>%
  mutate(skill = factor(skill),
         cohort = cohort*5+15,
         mod = case_when(
           mod == 'base' ~ "Base",
           mod == 'FC' ~ "Form. Cont.",
           mod == 'AS' ~ "Form. Preretr.",
         ),
         mod = factor(mod)) %>%
  mutate(value_pct = 100*(value/value[mod == 'Base'] - 1))

## plots in levels

plot_prodind_ss_w_lvl <- prodind_ss %>%
  filter(vars == "w") %>%
  ggplot(aes(x = cohort, y = value, colour = mod, linetype = skill)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_prodind_ss_c_lvl <- prodind_ss %>%
  filter(vars == "c" & variable != 'c') %>%
  ggplot(aes(x = cohort, y = value, colour = mod, linetype = skill)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_prodind_ss_k_lvl <- prodind_ss %>%
  filter(vars == "k") %>%
  ggplot(aes(x = cohort, y = value, colour = mod, linetype = skill)) +
  geom_hline(yintercept=0)+
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_prodind_ss_welf_lvl <- prodind_ss %>%
  filter(vars == "welf") %>%
  ggplot(aes(x = cohort, y = value, colour = mod, linetype = skill)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_prodind_ss_lab_lvl <- prodind_ss %>%
  filter(vars == "lab") %>%
  ggplot(aes(x = cohort, y = value, colour = mod, linetype = skill)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_prodind_ss_pen_lvl <- prodind_ss %>%
  filter(vars == "pen") %>%
  ggplot(aes(x = cohort, y = value, colour = mod, linetype = skill)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")


ggsave(filename = "tr2_prodind_ss_lvl_wage.eps",
       plot = plot_prodind_ss_w_lvl,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr2_prodind_ss_lvl_conso.eps",
       plot = plot_prodind_ss_c_lvl,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr2_prodind_ss_lvl_k.eps",
       plot = plot_prodind_ss_k_lvl,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr2_prodind_ss_lvl_welf.eps",
       plot = plot_prodind_ss_welf_lvl,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr2_prodind_ss_lvl_labsupp.eps",
       plot = plot_prodind_ss_lab_lvl,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr2_prodind_ss_lvl_pen.eps",
       plot = plot_prodind_ss_pen_lvl,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)


## plots in pct

plot_prodind_ss_w_pct <- prodind_ss %>%
  filter(vars == "w" & mod != 'Base') %>%
  ggplot(aes(x = cohort, y = value_pct, colour = mod, linetype = skill)) +
  geom_hline(yintercept=0) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_prodind_ss_c_pct <- prodind_ss %>%
  filter(vars == "c" & variable != 'c' & mod != 'Base') %>%
  ggplot(aes(x = cohort, y = value_pct, colour = mod, linetype = skill)) +
  geom_hline(yintercept=0) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_prodind_ss_k_pct <- prodind_ss %>%
  filter(vars == "k" & mod != 'Base') %>%
  ggplot(aes(x = cohort, y = value_pct, colour = mod, linetype = skill)) +
  geom_hline(yintercept=0) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_prodind_ss_welf_pct <- prodind_ss %>%
  filter(vars == "welf" & mod != 'Base') %>%
  ggplot(aes(x = cohort, y = value_pct, colour = mod, linetype = skill)) +
  geom_hline(yintercept=0) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_prodind_ss_lab_pct <- prodind_ss %>%
  filter(vars == "lab" & mod != 'Base') %>%
  ggplot(aes(x = cohort, y = value_pct, colour = mod, linetype = skill)) +
  geom_hline(yintercept=0) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_prodind_ss_pen_pct <- prodind_ss %>%
  filter(vars == "pen" & mod != 'Base') %>%
  ggplot(aes(x = cohort, y = value_pct, colour = mod, linetype = skill)) +
  geom_hline(yintercept=0) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")


ggsave(filename = "tr2_prodind_ss_pct_wage.eps",
       plot = plot_prodind_ss_w_pct,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr2_prodind_ss_pct_conso.eps",
       plot = plot_prodind_ss_c_pct,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr2_prodind_ss_pct_k.eps",
       plot = plot_prodind_ss_k_pct,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr2_prodind_ss_pct_welf.eps",
       plot = plot_prodind_ss_welf_pct,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr2_prodind_ss_pct_labsupp.eps",
       plot = plot_prodind_ss_lab_pct,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr2_prodind_ss_pct_pen.eps",
       plot = plot_prodind_ss_pen_pct,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

# tables in lvl

prodind_ss_lvl_macro <- prodind_ss %>%
  filter(variable %in% c('y', 'nbar', 'c', 'L_Q', 'L_NQ', 'kbar', 'H', 'I', 'R', 'rd')) %>%
  select(-variable, -cohort, -value_pct) %>%
  pivot_wider(names_from = mod,
              values_from = value)

prodind_ss_lvl_macro %>%
  xtable() %>%
  print(type = 'latex',
        file = '../output/tr2_prodind_ss_macro_lvl.tex')

prodind_ss_lvl_fiscal <- prodind_ss %>%
  filter(grepl('ratio', vars) |
           grepl('tau', vars) |
           grepl('T', vars) |
           vars %in% c("g", "D", "penbase", "penind", "retire", "Deped")) %>%
  select(mod, vars, skill, value, -cohort, -variable) %>%
  pivot_wider(names_from = mod,
              values_from = value)

prodind_ss_lvl_fiscal %>%
  xtable() %>%
  print(type = 'latex',
        file = '../output/tr2_prodind_ss_fiscal_lvl.tex')


## tables in pct
prodind_ss_pct_macro <- prodind_ss %>%
  filter(variable %in% c('y', 'nbar', 'c', 'L_Q', 'L_NQ', 'kbar', 'H', 'I', 'R', 'rd')) %>%
  select(-variable, -cohort, -value) %>%
  pivot_wider(names_from = mod,
              values_from = value_pct) %>%
  select(-Base)

prodind_ss_pct_macro %>%
  xtable() %>%
  print(type = 'latex',
        file = '../output/tr2_prodind_ss_macro_pct.tex')

prodind_ss_pct_fiscal <- prodind_ss %>%
  filter(grepl('ratio', vars) |
           grepl('tau', vars) |
           grepl('T', vars) |
           vars %in% c("g", "D", "penbase", "penind", "retire", "Deped")) %>%
  select(mod, vars, skill, value_pct, -cohort, -variable) %>%
  pivot_wider(names_from = mod,
              values_from = value_pct )%>%
  select(-Base)

prodind_ss_pct_fiscal %>%
  xtable() %>%
  print(type = 'latex',
        file = '../output/tr2_prodind_ss_fiscal_pct.tex')



# TR2: dynamic effects of individual productivities -----------------------


prodind_dyn <- read_csv('../output/tr2_skill_dyn.csv',
                        show_col_types = F) %>%
  separate(variable,
           into = c("vars", "skill", "cohort"),
           sep = '_',
           remove = F,
           convert = T) %>%
  mutate(skill = factor(skill),
         cohort = cohort*5+15,
         period = period*5+1900,
         mod = case_when(
           mod == 'base' ~ "Base",
           mod == 'fc' ~ "Form. Cont.",
           mod == 'as' ~ "Form. Preretr.",
         ),
         mod = factor(mod)) %>%
  mutate(value_pct = 100*(value/value[mod == 'Base'] - 1))


## plots in level

plot_prodind_dyn_y_lvl <- prodind_dyn %>%
  filter(vars == 'y') %>%
  ggplot(aes(x = period, y = value, colour = mod)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_prodind_dyn_nbar_lvl <- prodind_dyn %>%
  filter(vars == 'nbar') %>%
  ggplot(aes(x = period, y = value, colour = mod)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_prodind_dyn_lab_lvl <- prodind_dyn %>%
  filter(vars == 'L') %>%
  ggplot(aes(x = period,
             y = value,
             colour = mod,
             group = interaction(mod, skill),
             linetype = skill)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")



ggsave(filename = "tr2_prodind_dyn_lvl_y.eps",
       plot = plot_prodind_dyn_y_lvl,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr2_prodind_dyn_lvl_nbar.eps",
       plot = plot_prodind_dyn_nbar_lvl,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr2_prodind_dyn_lvl_labour.eps",
       plot = plot_prodind_dyn_lab_lvl,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)


## plots in pct

plot_prodind_dyn_y_pct <- prodind_dyn %>%
  filter(vars == 'y' & mod != 'Base') %>%
  ggplot(aes(x = period, y = value_pct, colour = mod)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_prodind_dyn_nbar_pct <- prodind_dyn %>%
  filter(vars == 'nbar' & mod != 'Base') %>%
  ggplot(aes(x = period, y = value_pct, colour = mod)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_prodind_dyn_lab_pct <- prodind_dyn %>%
  filter(vars == 'L' & mod != 'Base') %>%
  ggplot(aes(x = period,
             y = value_pct,
             colour = mod,
             group = interaction(mod, skill),
             linetype = skill)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")



ggsave(filename = "tr2_prodind_dyn_pct_y.eps",
       plot = plot_prodind_dyn_y_pct,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr2_prodind_dyn_pct_nbar.eps",
       plot = plot_prodind_dyn_nbar_pct,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr2_prodind_dyn_pct_labour.eps",
       plot = plot_prodind_dyn_lab_pct,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

# tables in pct

prodind_dyn %>%
  filter(vars %in% c('nbar', 'y', 'H', 'L', 'Ptot') &
           period %in% yrs &
           mod != 'Base') %>%
  select(-cohort, -value, -variable) %>%
  pivot_wider(names_from = 'period',
              values_from = 'value_pct') %>%
  arrange(mod, desc(vars), skill) %>%
  xtable() %>%
  print(type = 'latex', file = '../output/tr2_prodind_dyn.tex')



# TR2: skill incentive +5% ------------------------------------------------

skill_dyn <- read_csv('../output/tr2_skilledshock_dyn.csv',
                      show_col_types = F) %>%
  separate(variable,
           into = c("vars", "skill", "cohort"),
           sep = '_',
           remove = F,
           convert = T) %>%
  mutate(skill = factor(skill),
         cohort = cohort*5+15,
         period = period*5+1900,
         mod = case_when(
           mod == 'base' ~ "Base",
           mod == 'shock' ~ "Choc +5%",
         ),
         mod = factor(mod)) %>%
  mutate(value_pct = 100*(value/value[mod == 'Base'] - 1))

# plots in levels

plot_skill_dyn_y_lvl <- skill_dyn %>%
  filter(vars == 'y') %>%
  ggplot(aes(x = period, y = value, colour = mod)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_skill_dyn_L_lvl <- skill_dyn %>%
  filter(vars == 'L') %>%
  ggplot(aes(x = period, y = value, colour = mod, linetype = skill)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_skill_dyn_nbar_lvl <- skill_dyn %>%
  filter(vars == 'nbar') %>%
  ggplot(aes(x = period, y = value, colour = mod)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")

plot_skill_dyn_ptot_lvl <- skill_dyn %>%
  filter(vars == 'Ptot') %>%
  ggplot(aes(x = period, y = value, colour = mod)) +
  geom_line() + theme_minimal() + xlab('') + ylab('') +
  theme(legend.title = element_blank(), legend.position="bottom")


ggsave(filename = "tr2_skill_dyn_lvl_pib.eps",
       plot = plot_skill_dyn_y_lvl,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)


ggsave(filename = "tr2_skill_dyn_lvl_labour.eps",
       plot = plot_skill_dyn_L_lvl,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)


ggsave(filename = "tr2_skill_dyn_lvl_nbar.eps",
       plot = plot_skill_dyn_nbar_lvl,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)

ggsave(filename = "tr2_skill_dyn_lvl_poptot.eps",
       plot = plot_skill_dyn_ptot_lvl,
       device = 'eps',
       path = "../plots/",
       units = 'in',
       dpi = 'retina',
       width = 8,
       height = 4.5)


# tables in pct

skill_dyn %>%
  filter(vars %in% c('nbar', 'y', 'H', 'L', 'Ptot') &
           period %in% yrs &
           mod != 'Base') %>%
  select(-cohort, -value, -variable) %>%
  pivot_wider(names_from = 'period',
              values_from = 'value_pct') %>%
  arrange(mod, desc(vars), skill) %>%
  select(-mod) %>%
  xtable() %>%
  print(type = 'latex', file = '../output/tr2_skill_dyn.tex')