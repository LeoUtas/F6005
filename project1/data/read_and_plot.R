rm(list = ls())

setwd("D:/OneDrive - University of Tasmania/HOANG.Ng84/Education/Hoang.MUN20/Required courses/Fish 6005/labs/F6005/Project1_mess3/data")

library(ggplot2)
library(ggridges)
library(ggpubr)

vec_func <- function(x) {
  vec_data <- tidyr::gather(x, "Year", "Catch", -Length)
  return(vec_data)
}

ns <- 4

survey <- list()
for (i in 1:ns) {
  CL <- openxlsx::read.xlsx("data.xlsx", sheet = 2)
  survey[[i]] <- openxlsx::read.xlsx("data.xlsx", sheet = i + 4)
}

# ********* landings ********************************;

landings <- read.table("landings.txt", header = TRUE)
landings$Catch <- landings$Catch / 1000
landings$TAC <- landings$TAC / 1000

p <- ggplot(landings) +
  geom_line(aes(x = Year, y = Catch), size = 1) +
  geom_point(aes(x = Year, y = TAC), color = "red") +
  labs(x = "Year", y = "Landings (Kt)") +
  theme_minimal()

gname <- "landings_ts.jpeg"
jpeg(file = gname, width = 4, height = 4, units = "in", res = 300)
print(p)
dev.off()

################## Catch @ Length ####################################################;

CL_vec <- vec_func(CL)
CL_vec$Year <- as.numeric(CL_vec$Year)
CL_vec$symbol <- 19
CL_vec$symbol[CL_vec$Catch == 0] <- 8
CL_vec$size <- sqrt(CL_vec$Catch)
CL_vec$size[CL_vec$Catch == 0] <- 5

C_bubbles <- ggplot(CL_vec, aes(x = Year, y = Length)) +
  geom_point(aes(size = size, shape = factor(symbol)), alpha = 0.6, show.legend = FALSE) +
  scale_shape_manual(values = c(8, 19)) +
  scale_size_continuous(range = c(0, 8)) + 
  theme_minimal()
  

C_ridges <- ggplot(CL_vec) +
  theme_ridges(grid = T, center_axis_labels = T, font_size = 10) +
  geom_density_ridges(aes(x = Length, y = Year, height = Catch, fill = as.factor(Year)),
    stat = "identity", scale = 3, alpha = .8, show.legend = F
  ) +
  scale_y_continuous(breaks = 1991:2019, expand = c(0, 0)) +
  scale_x_continuous(limits = c(10, 45), breaks = c(10, 15, 20, 25, 30, 35, 40, 45), expand = c(0, 0)) +
  xlab("Length (cm)") +
  ylab("Year") +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 10, r = 20, b = 0, l = 0)),
    plot.margin = unit(c(0.25, 1, 0.5, 0.25), "cm"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.75)
  )

jpeg(file = "catch.jpeg", width = 16, height = 8, units = "in", res = 300)
p = ggarrange(C_bubbles, C_ridges)
print(p)
dev.off()


################## FALL RV @ Length ####################################################;
survey_vec <- lapply(survey, vec_func)

S_bubbles <- list()
S_ridges <- list()
plist = list()
id <- list(
  "Fall_RV",
  "Spring_RV",
  "Spanish_summer",
  "Spanish_3N"
)

for (i in 1:ns) {
  survey_vec[[i]]$Year <- as.numeric(survey_vec[[i]]$Year)
  survey_vec[[i]]$symbol <- 19
  survey_vec[[i]]$symbol[survey_vec[[i]]$Catch == 0] <- 8
  survey_vec[[i]]$size <- sqrt(survey_vec[[i]]$Catch)
  survey_vec[[i]]$size[survey_vec[[i]]$Catch == 0] <- 25

  S_bubbles[[i]] <- ggplot(survey_vec[[i]], aes(x = Year, y = Length)) +
    geom_point(aes(size = size, shape = factor(symbol)), alpha = 0.6, show.legend = FALSE) +
    scale_shape_manual(values = c(8, 19)) +
    scale_size_continuous(range = c(0, 10)) +
    theme_minimal()

  S_ridges[[i]] <- ggplot(survey_vec[[i]]) +
    theme_ridges(grid = T, center_axis_labels = T, font_size = 10) +
    geom_density_ridges(aes(x = Length, y = Year, height = Catch, fill = as.factor(Year)),
      stat = "identity", scale = 5, alpha = .8, show.legend = F
    ) +
    xlab("Length (mm)") +
    ylab("Year") +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
      axis.title.x = element_text(margin = margin(t = 10, r = 20, b = 0, l = 0)),
      plot.margin = unit(c(0.25, 1, 0.5, 0.25), "cm"),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.75)
    )
  
  jpeg(filename = paste(id[[i]], "survey.jpeg", sep = "_"), width = 16, height = 8, units = "in", res = 300)
  plist[[i]] = ggarrange(S_bubbles[[i]], S_ridges[[i]])
  print(plist[[i]])
  dev.off()
  }

############### maturity gives #######

mat <- read.table("maturity.txt", header = TRUE)

gname <- "maturity_ogives.jpeg"
p <- ggplot(mat) +
  geom_line(aes(x = Length, y = P1, color = "black"), size = 1) +
  geom_line(aes(x = Length, y = P2, color = "red"), size = 1) +
  labs(x = "Length", y = "Female Proportion Mature") +
  scale_colour_manual(
    name = "Time Periods",
    values = c("black" = "black", "red" = "red"), labels = c("1991-1995", "1996+")
  ) +
  theme_minimal()

jpeg(file = gname, width = 10, height = 8, units = "in", res = 300)
print(p)
dev.off()

############### Weight ~ length #######

year <- 1990:2019
len <- 1:65
pdat <- expand.grid(year = year, length = len)

LW_parms <- read.table("LW_parms.txt", header = TRUE)

pdat$a <- rep(LW_parms$a, times = length(len))
pdat$b <- rep(LW_parms$b, times = length(len))

pdat$weight <- pdat$a * (pdat$length**pdat$b) / 1000

library(viridis)
library(hrbrthemes)

personal_theme <- theme(
  
  plot.title = element_text(hjust = .5),
  axis.title.x = element_text(hjust = .5),
  axis.title.y = element_text(hjust = .5)
)

p <- ggplot(pdat, aes(x = length, y = weight, color = as.factor(year))) +
  geom_line() +
  scale_color_viridis(discrete = TRUE) +
  ggtitle("Weight ~ Length") +
  theme_ipsum() +
  labs(x = "Length(cm)", y = "Weight (Kg)", color = "Year") + 
  theme_minimal()

gname <- "LW.jpeg"
jpeg(file = gname, width = 6, height = 4, units = "in", res = 300)
p + personal_theme
dev.off()

wt_mat <- matrix(pdat$weight,
  nrow = length(len),
  ncol = length(year),
  byrow = TRUE
)
rownames(wt_mat) <- len
colnames(wt_mat) <- year

wt_vec <- pdat[, names(pdat) %in% c("year", "length", "weight")]

save(
  landings,
  CL,
  CL_vec,
  survey,
  survey_vec,
  mat,
  wt_mat,
  wt_vec,
  LW_parms,
  file = "3LN_redfish.RData"
)
