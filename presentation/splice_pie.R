#!/usr/bin/env Rscript
library(ggplot2)
library(ggpubr)

input <- file('stdin', 'r')
ase_number <- readLines(input, n = 1)
cse_number <- readLines(input, n = 1)
ase_total <- readLines(input, n = 1)
cse_total <- readLines(input, n = 1)

ase_cse <-
  matrix(c(
    as.numeric(cse_total),
    as.numeric(ase_total),
    as.numeric(cse_number),
    as.numeric(ase_number)
  ), nrow = 2)
test <- chisq.test(ase_cse)

args <- commandArgs(T)
df <- data.frame(group = c("ASE", "CSE"),
                 value = c(as.numeric(ase_number), as.numeric(cse_number)))
labs <-
  paste0(df$group,
         "\n",
         df$value,
         "(",
         paste(round(df$value / sum(df$value), 4) * 100, "%", sep = ''),
         ")")
p <-
  ggpie(
    df,
    "value",
    label = labs,
    lab.font = c("bold", "#000000"),
    fill = "group",
    font.family = "ArialMT",
    palette = c("#7FA8EA", "#CC6666"),
    ggtheme = theme_transparent(base_size = 5),
    color = "#000000",
    size = 0.11672,
    legend = "NA"
  )

result_save_path <- args[1]
ggsave(
  result_save_path,
  plot = p,
  device = "pdf",
  width = .45,
  height = .45,
  dpi = "retina"
)

a <- prettyNum(as.numeric(cse_total),
               big.mark = ",",
               scientific = F)
b <- prettyNum(as.numeric(ase_total),
               big.mark = ",",
               scientific = F)
c <- prettyNum(as.numeric(cse_number),
               big.mark = ",",
               scientific = F)
d <- prettyNum(as.numeric(ase_number),
               big.mark = ",",
               scientific = F)
ase_densi <- format(
  as.numeric(ase_number) / as.numeric(ase_total),
  digits = 2,
  scientific = T
)
cse_densi <- format(
  as.numeric(cse_number) / as.numeric(cse_total),
  digits = 2,
  scientific = T
)
pvalue <- format(test[["p.value"]], digits = 2, scientific = T)
pdf(args[2], width = 270 / 127, height = 135 / 127)
par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(
  NA,
  xlim = c(0, 270 / 127),
  ylim = c(0, 135 / 127),
  bty = 'n',
  xaxt = 'n',
  yaxt = 'n',
  xlab = '',
  ylab = ''
)
text(19.5 * 5 / 127,
     18.075 * 5 / 127,
     "Total",
     pos = 2,
     cex = 5 / 12)
text(26 * 5 / 127,
     18.075 * 5 / 127,
     "Nm",
     pos = 2,
     cex = 5 / 12)
text(8.6 * 5 / 127,
     14.446 * 5 / 127,
     "CSE",
     pos = 2,
     cex = 5 / 12)
text(8.6 * 5 / 127,
     11.386 * 5 / 127,
     "ASE",
     pos = 2,
     cex = 5 / 12)
text(8.6 * 5 / 127,
     8.326 * 5 / 127,
     "P value",
     pos = 2,
     cex = 5 / 12)
text(19.5 * 5 / 127, 14.446 * 5 / 127, a, pos = 2, cex = 5 / 12)
text(19.5 * 5 / 127, 11.386 * 5 / 127, b, pos = 2, cex = 5 / 12)
text(26 * 5 / 127, 14.446 * 5 / 127, c, pos = 2, cex = 5 / 12)
text(26 * 5 / 127, 11.386 * 5 / 127, d, pos = 2, cex = 5 / 12)
text(35 * 5 / 127,
     14.446 * 5 / 127,
     cse_densi,
     pos = 2,
     cex = 5 / 12)
text(35 * 5 / 127,
     11.386 * 5 / 127,
     ase_densi,
     pos = 2,
     cex = 5 / 12)
text(19.5 * 5 / 127,
     8.326 * 5 / 127,
     pvalue,
     pos = 2,
     cex = 5 / 12)
lines(x = c(6.8 * 5 / 127, 6.8 * 5 / 127),
      y = c(19.5 * 5 / 127, 10.5 * 5 / 127))
lines(x = c(2 * 5 / 127, 23 * 5 / 127),
      y = c(16.657 * 5 / 127, 16.657 * 5 / 127))