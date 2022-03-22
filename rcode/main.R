### General script for stuff with CDC model simulations

library(tidyverse)
library(R.matlab)
library(matlabr)
library(xtable)
options(xtable.include.rownames = F)



# historical life profiles ------------------------------------------------

pha <- read_csv(file = '../noyau_cdc_PHA/output/pha_baseline.csv')

pha_ss <- pha[1, ] %>%
  pivot_longer(cols = names(.),
               names_to = "vars",
               values_to = 'ss') %>%
  separate(vars,
           into = c('variable', 'qual', 'cohort'),
           sep = '_', remove = F)

pha_sims <- pha[c(1, 241:280), ] %>%
  mutate(period = c(0, 1900 + (1:(n()-1))*5)) %>%
  pivot_longer(cols = -period,
               names_to = "vars",
               values_to = "value") %>%
  separate(vars,
           into = c('variable', 'qual', 'cohort'),
           sep = '_', remove = F) %>%
  mutate(qual = factor(qual, levels = c("Q", "NQ")),
         cohort = as.numeric(cohort))

plot_beta_pha <- pha_sims %>%
  filter(variable %in% c("beta") &
           period %in% c(1905, 1950, 2000, 2050, 0)) %>%
  arrange(period, cohort, variable) %>%
  ggplot(aes(x = cohort, y = value, group = period, colour = as.factor(period))) +
  geom_line() +
  facet_wrap(.~ qual)

plot_h_pha <- pha_sims %>%
  filter(variable %in% c("h") &
           period %in% c(1905, 1950, 2000, 2050, 0)) %>%
  arrange(period, cohort, variable) %>%
  ggplot(aes(x = cohort, y = value, group = period, colour = as.factor(period))) +
  geom_line() +
  facet_wrap(.~ qual, scales = 'free')

plot_p_pha <- pha_sims %>%
  filter(variable %in% c("P") &
           period %in% c(1905, 1950, 2000, 2050, 0)) %>%
  arrange(period, cohort, variable) %>%
  ggplot(aes(x = cohort, y = value, group = period, colour = as.factor(period))) +
  geom_line() +
  facet_wrap(.~ qual)



# Steady state comparisons full CDC & pension reforms --------------------------

# find matlab and run the script generating the files
get_matlab()
run_matlab_script("../cdc_maquette_ss_pensprodind/cdc_full_prodindpens.m",
                  verbose = F)

pens_ss <- read_csv(file = "../cdc_maquette_ss_pensprodind/output/pens_ss.csv",
                    show_col_types = F) %>%
  pivot_wider(names_from = 'model', values_from = "value") %>%
  separate(col = "variable",
           into = c("vars", "skill", "cohort"),
           sep = "_",
           remove = F) %>%
  mutate(skill = factor(skill, c("NQ", "Q")),
         cohort = as.numeric(cohort)*5 +15)

pens_vars <- pens_ss %>%
  filter(vars %in% c("welf", "w",
                     "P", "k", 'lab',
                     "Tw", "Tc", "Tk",
                     "Twratio", "Tcratio", "Tkratio",
                     "tauw", "tauc", "tauk", "tauf",
                     "Gratio", "g",
                     "Def", "Defratio",
                     "D", "Dratio",
                     "retire", "Penratio",
                     "pen", "penind", "penbase",
                     "Pret",
                     "Edratio", "Deped",
                     "y", "nbar", "I", "c")) %>%
  mutate(across(.cols = c("55", "65", "70"),
                .names = "{.col}_pct",
                .fns = ~ (.x-base)*100/base))

pens_fiscal <- pens_vars %>%
  select(-c(variable,
            base, `55`, `65`, `70`,
            skill, cohort)) %>%
  filter(grepl('ratio', vars) |
           grepl('tau', vars) |
           grepl('T', vars) |
           vars %in% c("g", "D", "penbase", "penind", "retire", "Deped"))

pens_econ <- pens_vars %>%
  select(vars, contains("_pct")) %>%
  filter(vars %in% c("y", "nbar", "I"))

pens_econ %>%
  xtable() %>%
  print(type = 'latex',
        file = "../output/tab_ss_pens_econ.tex",
        floating = F,
        latex.environment = NULL)


pens_fiscal %>%
  xtable() %>%
  print(type = 'latex',
        file = "../output/tab_ss_pens_fiscal.tex",
        floating = F,
        latex.environment = NULL)

plot_pens_welf <- pens_vars %>%
  filter(vars == 'welf') %>%
  select(-c("base", `55`, `65`, `70`)) %>%
  pivot_longer(cols = contains("_pct"), names_to = "model", values_to = 'value') %>%
  ggplot(aes(x = cohort, y = value, colour = model)) +
  geom_line() +
  facet_wrap(.~skill, nrow = 2) +
  ylab('Pourcentage par rapport au depart à 60 ans') +
  xlab("Age") +
  theme_bw() +
  ggtitle('Bienetre')

plot_pens_conso <- pens_vars %>%
  filter(vars == 'c' & variable != 'c') %>%
  select(-c("base", `55`, `65`, `70`)) %>%
  pivot_longer(cols = contains("_pct"), names_to = "model", values_to = 'value') %>%
  ggplot(aes(x = cohort, y = value, colour = model)) +
  geom_line() +
  facet_wrap(.~skill, nrow = 2) +
  ylab('Pourcentage par rapport au depart à 60 ans') +
  xlab("Age") +
  theme_bw() +
  ggtitle('Consommation')


plot_pens_montret <- pens_vars %>%
  filter(vars == 'pen') %>%
  select(-c("base", `55`, `65`, `70`)) %>%
  pivot_longer(cols = contains("_pct"), names_to = "model", values_to = 'value') %>%
  ggplot(aes(x = cohort, y = value, colour = model)) +
  geom_line() +
  facet_wrap(.~skill, nrow = 2) +
  ylab('Pourcentage par rapport au depart à 60 ans') +
  xlab("Age") +
  theme_bw() +
  ggtitle('Montant des Retraites')

ggsave(plot = plot_pens_welf,
       path = "../plots/",
       filename = "pens_welf.eps",
       device = 'eps',
       dpi = 'retina')

ggsave(plot = plot_pens_conso,
       path = "../plots/",
       filename = "pens_conso.eps",
       device = 'eps',
       dpi = 'retina')

ggsave(plot = plot_pens_montret,
       path = "../plots/",
       filename = "pens_montret.eps",
       device = 'eps',
       dpi = 'retina')


plot_pens_pop <- pens_vars %>%
  filter(vars == 'k') %>%
  select(-c("base", `55`, `65`, `70`)) %>%
  pivot_longer(cols = contains("_pct"), names_to = "model", values_to = 'value') %>%
  ggplot(aes(x = cohort, y = value, colour = skill, linetype = model)) +
  geom_line()



# SS comparison: full CDC, individual productivities ----------------------

prodind <- read_csv()
