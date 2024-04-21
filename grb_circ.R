grb = read.csv("grb_circular_dataset.csv")
head(grb)
library(circular)
library(CircStats)

hist(grb$T90.s)
hist(log(grb$T50.s))
summary(log(grb$T90))

library(circular)
library(Directional)
library(CircStats)

watson.test(grb$GLON.deg, alpha = 0.01, dist = "vonmises")
watson.test(grb$GLAT.deg, alpha = 0.01, dist = "vonmises")

vmf_density_grid = 
  function(u, ngrid = 100) {
    # Translate to (0,180) and (0,360)
    u[,1] <- u[,1] + 90
    u[,2] <- u[,2] + 180
    res <- vmf.kerncontour(u, thumb = "none", den.ret = T, full = T,
                           ngrid = ngrid)
    
    # Translate back to (-90, 90) and (-180, 180) and create a grid of
    # coordinates
    ret <- expand.grid(Lat = res$lat - 90, Long = res$long - 180)
    ret$Density <- c(res$den)
    ret
  }

grb_c = cbind(grb$GLAT.deg , grb$GLON.deg)
set.seed(2022)
EvMFs <- 
  function(K){
    movMF(grb_c, k = K, control= list(nruns = 20))
  }

Esd = lapply(1:10, EvMFs)
gt = sapply(Esd, BIC)
gt

Esd

library(ggplot2)
#Plot histogram of T90 values for short and long GRBs
ggplot(grb, aes(x = log(grb$T50.s), fill = type)) +
  geom_histogram(binwidth = 1, color = "black", alpha = 0.7) +
  labs(x = "log transformed T50 (seconds)", y = "Frequency", title = "Distribution of T50 for Short and Long GRBs") +
  scale_fill_manual(values = c("Short" = "blue", "Long" = "red")) +
  theme_minimal()

ggplot(grb, aes(x = GLON.deg, y = GLAT.deg, color = type)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Long" = "red", "Short" = "blue")) +
  labs(x = "Galactic Longitude", y = "Galactic Latitude", title = "Sky Map of GRBs") +
  theme_minimal() +
  coord_map("aitoff")

?movMF::dmovMF()


long_grbs <- grb[grb$type == "Long", ]
short_grbs <- grb[grb$type == "Short", ]

#long_grbs

#grb

watson.test(long_grbs$GLON.deg, alpha = 0.01, dist = "vonmises")
watson.test(long_grbs$GLAT.deg, alpha = 0.01, dist = "vonmises")

grb_L = cbind(long_grbs$GLAT.deg, long_grbs$GLON.deg)
set.seed(2023)
EvMFs <- 
  function(K){
    movMF(grb_L, k = K, control= list(nruns = 20))
  }
Esd = lapply(1:10, EvMFs)
gt = sapply(Esd, BIC)
gt



watson.test(short_grbs$GLON.deg, alpha = 0.01, dist = "vonmises")
watson.test(short_grbs$GLAT.deg, alpha = 0.01, dist = "vonmises")


grb_S = cbind(short_grbs$GLAT.deg, short_grbs$GLON.deg)
set.seed(2024)
EvMFs <- 
  function(K){
    movMF(grb_S, k = K, control= list(nruns = 20))
  }
Esd = lapply(1:10, EvMFs)
gt = sapply(Esd, BIC)
gt

GRB.densities <- vmf_density_grid(grb[,c("GLAT.deg",
                                                    "GLON.deg")],
                                       ngrid = 300);

g.batse <- ggplot(grb, aes(x = GLON.deg, y = GLAT.deg)) +
  geom_point(size = 1, color = "black") +
  geom_density_2d(data = grb,
                  aes(x = GLON.deg, y = GLAT.deg),
                  color = "green", alpha = 1) +
  geom_contour(data = GRB.densities, aes(x=Long, y=Lat, z=Density),
               color = "blue") +
  theme_minimal() +
  coord_map("aitoff")

g.batse

Long_GRB.densities <- vmf_density_grid(long_grbs[,c("GLAT.deg",
                                                   "GLON.deg")],
                                       ngrid = 300);



g.long <- ggplot(long_grbs, aes(x = GLON.deg, y = GLAT.deg)) +
  geom_point(size = 1, color = "black") +
  geom_density_2d(data = long_grbs,
                  aes(x = GLON.deg, y = GLAT.deg),
                  color = "green", alpha = 1) +
  geom_contour(data = Long_GRB.densities, aes(x=Long, y=Lat, z=Density),
               color = "blue") +
  theme_minimal() +
  coord_map("aitoff")
g.long




Short_GRB.densities <- vmf_density_grid(short_grbs[,c("GLAT.deg",
                                                    "GLON.deg")],
                                       ngrid = 300);

g.short <- ggplot(short_grbs, aes(x = GLON.deg, y = GLAT.deg)) +
  geom_point(size = 1, color = "black") +
  geom_density_2d(data = short_grbs,
                  aes(x = GLON.deg, y = GLAT.deg),
                  color = "green", alpha = 1) +
  geom_contour(data = Short_GRB.densities, aes(x=Long, y=Lat, z=Density),
               color = "blue") +
  theme_minimal() +
  coord_map("aitoff")
g.short


feri = read.csv("fermi_location.csv")
head(feri)
watson.test(feri$RAJ2000, alpha = 0.01, dist = "vonmises")
watson.test(feri$DEJ2000, alpha = 0.01, dist = "vonmises")

feri_mat = cbind(feri$DEJ2000, feri$RAJ2000)
fermi.densities <- vmf_density_grid(feri[,c("DEJ2000",
                                                    "RAJ2000")], ngrid = 300)
EvMFs <- 
  function(K){
    movMF(feri_mat, k = K, control= list(nruns = 20))
  }
Esd = lapply(1:10, EvMFs)
gt = sapply(Esd, BIC)
gt


g.fermi <- ggplot(feri, aes(x = RAJ2000, y = DEJ2000)) +
  geom_point(size = 1, color = "black") +
  geom_density_2d(data = feri,
                  aes(x = RAJ2000, y = DEJ2000),
                  color = "green", alpha = 1) +
  geom_contour(data = fermi.densities, aes(x=Long, y=Lat, z=Density),
               color = "blue") +
  theme_minimal() +
  coord_map("aitoff")
g.fermi


head(feri)

fermi_time = read.csv("fermi_time.csv")
head(fermi_time)

fermi_data_time<- fermi_time %>%
  mutate(Type = ifelse(fermi_time$T90 > threshold, "Long", "Short"))

head(fermi_data_time)
ggplot(fermi_data_time, aes(x = log(fermi_time$T90), fill = Type)) +
  geom_histogram(binwidth = 1, color = "black", alpha = 0.7) +
  labs(x = "log transformed T90 (seconds)", y = "Frequency", title = "Distribution of T90 for Short and Long GRBs") +
  scale_fill_manual(values = c("Short" = "blue", "Long" = "red")) +
  theme_minimal()

ggplot(fermi_data_time, aes(x = log(fermi_time$T50), fill = Type)) +
  geom_histogram(binwidth = 1, color = "black", alpha = 0.7) +
  labs(x = "log transformed T50 (seconds)", y = "Frequency", title = "Distribution of T50 for Short and Long GRBs") +
  scale_fill_manual(values = c("Short" = "blue", "Long" = "red")) +
  theme_minimal()


# Merge datasets based on common_param
fermi_dataset <- merge(feri, fermi_data_time, by = "Fermi", all = TRUE)

# Display merged dataset
head(fermi_dataset)

fermi_dataset_f = fermi_dataset[,c(1,3,5,6,10,11,15,18)]
head(fermi_dataset_f)

fermi_data_final<- fermi_dataset_f %>%
  mutate(Type = ifelse(fermi_dataset_f$T90 > threshold, "Long", "Short"))

head(fermi_data_final)
fermi_data_final = na.omit(fermi_data_final)
dim(fermi_data_final)

long_Fgrb <- fermi_data_final[fermi_data_final$Type == "Long", ]
short_Fgrb <- fermi_data_final[fermi_data_final$Type == "Short", ]

watson.test(long_Fgrb$RAJ2000, alpha = 0.01, dist = "vonmises")
watson.test(long_Fgrb$DEJ2000, alpha = 0.01, dist = "vonmises")

ggplot(fermi_data_final, aes(x = RAJ2000, y = DEJ2000, color = Type)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Long" = "red", "Short" = "blue")) +
  labs(x = "Galactic Longitude", y = "Galactic Latitude", title = "Sky Map of GRBs") +
  theme_minimal() +
  coord_map("aitoff")

feri_mat_l = cbind(long_Fgrb$DEJ2000, long_Fgrb$RAJ2000)
L.fermi.densities <- vmf_density_grid(long_Fgrb[,c("DEJ2000",
                                            "RAJ2000")], ngrid = 300)
EvMFs <- 
  function(K){
    movMF(feri_mat_l, k = K, control= list(nruns = 20))
  }
Esd = lapply(1:10, EvMFs)
gt = sapply(Esd, BIC)
gt


g.fermi.l <- ggplot(long_Fgrb, aes(x = RAJ2000, y = DEJ2000)) +
  geom_point(size = 1, color = "black") +
  geom_density_2d(data = long_Fgrb,
                  aes(x = RAJ2000, y = DEJ2000),
                  color = "green", alpha = 1) +
  geom_contour(data = L.fermi.densities, aes(x=Long, y=Lat, z=Density),
               color = "blue") +
  theme_minimal() +
  coord_map("aitoff")
g.fermi.l

head(long_Fgrb)

watson.test(short_Fgrb$RAJ2000, alpha = 0.01, dist = "vonmises")
watson.test(short_Fgrb$DEJ2000, alpha = 0.01, dist = "vonmises")

feri_mat_s = cbind(short_Fgrb$DEJ2000, short_Fgrb$RAJ2000)
S.fermi.densities <- vmf_density_grid(short_Fgrb[,c("DEJ2000",
                                                   "RAJ2000")], ngrid = 300)
EvMFs <- 
  function(K){
    movMF(feri_mat_s, k = K, control= list(nruns = 20))
  }
Esd = lapply(1:10, EvMFs)
gt = sapply(Esd, BIC)
gt


g.fermi.s <- ggplot(short_Fgrb, aes(x = RAJ2000, y = DEJ2000)) +
  geom_point(size = 1, color = "black") +
  geom_density_2d(data = short_Fgrb,
                  aes(x = RAJ2000, y = DEJ2000),
                  color = "green", alpha = 1) +
  geom_contour(data = S.fermi.densities, aes(x=Long, y=Lat, z=Density),
               color = "blue") +
  theme_minimal() +
  coord_map("aitoff")
g.fermi.s

head(long_Fgrb)
