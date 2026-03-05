library(dplyr)
library(ggplot2)
library(tidyr)

# Read dataframe
m <- read.csv("mvroot_MDcat_CI.csv")
m$Method <- gsub("CASTLES-Pro", "CoalBL", m$Method)
m$Calibrations <- factor(m$Calibrations)
m$ADbin <- cut(
  m$ad,
  breaks=c(0,0.25,0.5,0.75,1),
  include.lowest=TRUE
)
m$outgroup = factor(grepl("outgroup.1", m$Condition))
m$ratevar = factor(sub(".genes.*","",sub("outgroup.*.species.","",m$Condition)))
outgroup.labs <- c("With outgroup","No outgroup")
names(outgroup.labs) <- c(TRUE, FALSE)
ratevar.labs <- c("Low","Medium", "High")
names(ratevar.labs) <- c("5","1.5", "0.15")
m$param <- as.numeric(m$mu)
m$param_lower <- as.numeric(m$mu_lower)
m$param_upper <- as.numeric(m$mu_upper)

coal <- m %>% filter(grepl("CoalBL", Method))
con  <- m %>% filter(grepl("ConBL", Method))
coal <- coal %>%
  select(ratevar, Calibrations, Replicate, clade_id,
         ad, param_coal = param)

con <- con %>%
  select(ratevar, Calibrations, Replicate, clade_id,
         ad, param_con = param, param_lower,param_upper)

merged <- inner_join(coal, con, by = c("ratevar","Calibrations", "Replicate","clade_id","ad"))
merged <- merged %>% mutate(inside = param_coal >= param_lower & param_coal <= param_upper)
merged <- merged %>%
  group_by(ratevar, Calibrations, Replicate) %>%
  mutate(
    param_decile = ntile(param_con, 10)
  ) %>% ungroup()

merged$ADbin <- cut(merged$ad, breaks=c(0,0.25,0.5,0.75,1), include.lowest=TRUE)
merged$param_decile <- factor(merged$param_decile)

merged$inside <- factor(merged$inside)

summary_df <- merged %>%
  group_by(Calibrations, ADbin, ratevar, inside) %>%
  summarise(n = n(), .groups="drop") %>%
  group_by(Calibrations, ADbin, ratevar) %>%
  mutate(prop = n/sum(n)) %>%
  ungroup()

avg_df <- summary_df %>%
  filter(inside == "TRUE") %>%
  group_by(ratevar, Calibrations, ADbin) %>%
  summarise(avg_prop = mean(prop), .groups="drop")

ggplot(summary_df,aes(x=param_decile, y=prop, fill=inside)) +
  geom_col(width=0.85) +
  geom_hline(data = avg_df, aes(yintercept = avg_prop), color = "black", linewidth = 1, linetype=2) +
  facet_grid(Calibrations ~ ADbin) +
  scale_y_continuous(limits=c(0,1), labels=scales::percent, name="Proportion of branches") +
  scale_fill_manual(
    values=c("FALSE"="#F2A45E","TRUE"="#ACE1AF"),
    name="CoalBL inside ConBL 95% CI") +
  scale_x_discrete(name="ConBL branch length decile") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=0), legend.position="bottom")
ggsave("CI_branches_AD_MVroot.pdf",
       width=8, height=5.8)

ggplot(summary_df,aes(x=param_decile, y=prop, fill=inside)) +
  geom_col(width=0.85) +
  geom_hline(data = avg_df, aes(yintercept = avg_prop), color = "black", linewidth = 1, linetype=2) +
  facet_grid(Calibrations ~ ADbin) +
  scale_y_continuous(limits=c(0,1), labels=scales::percent, name="Proportion of branches") +
  scale_fill_manual(
    values=c("FALSE"="#F2A45E","TRUE"="#ACE1AF"),
    name="CoalBL inside ConBL 95% CI") +
  scale_x_discrete(name="ConBL mutation rate decile") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=0), legend.position="bottom")
ggsave("CI_mu_rates_AD_MVroot.pdf",
       width=8, height=5.8)

ggplot(summary_df,aes(x=param_decile, y=prop, fill=inside)) +
  geom_col(width=0.85) +
  geom_hline(data = avg_df, aes(yintercept = avg_prop), color = "black", linewidth = 1, linetype=2) +
  facet_grid(Calibrations ~ ADbin) +
  scale_y_continuous(limits=c(0,1), labels=scales::percent, name="Proportion of nodes") +
  scale_fill_manual(
    values=c("FALSE"="#F2A45E","TRUE"="#ACE1AF"),
    name="CoalBL inside ConBL 95% CI") +
  scale_x_discrete(name="ConBL divergence time decile") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=0), legend.position="none")
ggsave("CI_node_ages_AD_MVroot.pdf",
       width=8, height=5)

ggplot(summary_df,aes(x=factor(Calibrations), y=prop, fill=inside)) +
  geom_col(width=0.8) +
  facet_grid(ratevar ~ ADbin,labeller = labeller(ratevar = ratevar.labs)) +
  scale_y_continuous(limits=c(0,1), labels=scales::percent, name="Proportion of branches") +
  scale_x_discrete(name="Number of calibrations") +
  scale_fill_manual(
    values=c("FALSE"="#F2A45E","TRUE"="#ACE1AF"),
    name="CoalBL inside ConBL 95% CI") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=0), legend.position="none")
ggsave("mvroot_branches_CI_summary.pdf",width=4, height=4)

ggplot(summary_df,aes(x=factor(Calibrations), y=prop, fill=inside)) +
  geom_col(width=0.8) +
  facet_grid(ratevar ~ ADbin,labeller = labeller(ratevar = ratevar.labs)) +
  scale_y_continuous(limits=c(0,1), labels=scales::percent, name="Proportion of node ages") +
  scale_x_discrete(name="Number of calibrations") +
  scale_fill_manual(
    values=c("FALSE"="#F2A45E","TRUE"="#ACE1AF"),
    name="CoalBL inside ConBL 95% CI") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=0), legend.position="none")
ggsave("mvroot_ages_CI_summary.pdf",width=4, height=4)

ggplot(summary_df,aes(x=factor(Calibrations), y=prop, fill=inside)) +
  geom_col(width=0.8) +
  facet_grid(ratevar ~ ADbin,labeller = labeller(ratevar = ratevar.labs)) +
  scale_y_continuous(limits=c(0,1), labels=scales::percent, name="Proportion of mutation rates") +
  scale_x_discrete(name="Number of calibrations") +
  scale_fill_manual(
    values=c("FALSE"="#F2A45E","TRUE"="#ACE1AF"),
    name="CoalBL inside ConBL 95% CI") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=0), legend.position="none")
ggsave("mvroot_mu_rates_CI_summary.pdf",width=4, height=4)





