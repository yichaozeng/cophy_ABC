setwd('/home/u29/yichaozeng/Desktop')

# specify which parameters to look at
cols <- 2:5

# plot the prior and posterior parameter distributions
para_inference <- read.csv(file = paste(prefix, 'cophy_ABC_statistics/', folder_id_slash, 'prior_posterior_real.csv', sep = ''))

para_inference$inference <- para_inference$bayes
para_inference$inference[para_inference$inference == 'posterior-before-regression'] <- 'posterior-ABC-rejection'
para_inference$inference[para_inference$inference == 'posterior-after-regression'] <- 'posterior-ABC-regression'

para_inference <- para_inference[para_inference$inference == 'prior' | para_inference$inference == 'posterior-ABC-rejection' | para_inference$inference == 'posterior-ABC-regression', ]
para_inference$inference <- factor(para_inference$inference, levels = c('prior', 'posterior-ABC-rejection', 'posterior-ABC-regression'))

# plotting
para_inference <- para_inference[order(para_inference$inference, decreasing = F), ]
para_inference <- para_inference[!is.na(para_inference$inference), ]

# convert dataset for plotting
para_inference_regression <- para_inference[para_inference$inference == 'posterior-ABC-regression',]
rate_dat <- data.frame(
  rate = c(para_inference_regression$lambda_H, para_inference_regression$lambda_S, para_inference_regression$lambda_C, para_inference_regression$exp_H) / 13.75,
  event = c(rep('Host speciation', nrow(para_inference_regression)), rep('Symbiont speciation', nrow(para_inference_regression)), rep('Copeciation', nrow(para_inference_regression)), rep('Host expansion speciation', nrow(para_inference_regression)) )
)
rate_dat$event <- factor(rate_dat$event, levels = c('Host speciation', 'Symbiont speciation', 'Copeciation', 'Host expansion speciation'))

library(GGally)
library(rlang)
library(grid)

rate_plot <-   ggplot(rate_dat, aes(rate)) +
  geom_density(aes(fill = event, color = event), linewidth = 0.5, alpha = 0, kernel = 'gaussian') +
  labs(x = "events per Myr") +
  theme_minimal() +  # Clean layout
  theme(
    legend.position.inside = c(0.8, 0.8),  # Position legend within each density plot panel
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key.size = unit(0.4, "lines"),
    text=element_text(family = "DejaVu Sans"),
    panel.border = element_rect(colour = "black", fill=NA)
  )


plot_name <- "ABC_SMC_ex_ind_real_1.png"
png(paste(plot_folder, plot_name, sep = ''), width = 6, height = 3, units = 'in', pointsize = 12, res = 300)

print(rate_plot)
# grid.newpage()
# grid.draw(pairs)
# vp = viewport(x=.7, y=.65, width=.35, height=.3) ## control legend position
# pushViewport(vp)
# grid.draw(legend)

dev.off()
