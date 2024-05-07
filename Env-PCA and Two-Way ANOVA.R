### Load library
library(tidyverse)
library(ggplot2)
library(vegan)

### Read all data from EPD website
setwd("C:/Users/")
e1 = read.csv("marine_water_quality_DM.csv",check.names = F)
e1 %>% head
e2 = read.csv("marine_water_quality_MM.csv",check.names = F)
e2 %>% head  
e3 = read.csv("marine_water_quality_NM.csv",check.names = F)
e4 = read.csv("marine_water_quality_PM.csv",check.names = F)
e5 = read.csv("marine_water_quality_SM.csv",check.names = F)
e6 = read.csv("marine_water_quality_TM.csv",check.names = F)
e7 = read.csv("marine_water_quality_VM.csv",check.names = F)
e8 = read.csv("marine_water_quality_WM.csv",check.names = F)
e9 = read.csv("marine_water_quality_EM.csv",check.names = F)
env = rbind(e1,e2,e3,e4,e5,e6,e7,e8,e9)

### Handling of "N/A", ">0.5" and "<0.01" values
env1<-
  env %>% as_tibble() %>% 
  filter(Depth == "Surface Water") %>% 
  mutate(Year = Dates %>% str_split("/") %>% sapply('[', 1)) %>% 
  mutate(Month = Dates %>% str_split("/") %>% sapply('[', 2)) %>% 
  mutate(Sample = paste(Station, Year, Month, sep = "_")) %>% 
  select(32,6:29) %>% 
  mutate(VSS = `Volatile Suspended Solids` %>% str_replace("<0.5", "0.25")) %>% 
  mutate(NH3 = `Unionised Ammonia` %>% str_replace("<0.001", "0.0005")) %>% 
  mutate(turbidity = Turbidity) %>% 
  mutate(Total_P = `Total Phosphorus` %>% str_replace("<0.02", "0.01")) %>% 
  mutate(Total_N = `Total Nitrogen`) %>%
  mutate(Total_KN = `Total Kjeldahl Nitrogen` %>% str_replace("<0.05", "0.025")) %>%
  mutate(TIN = `Total Inorganic Nitrogen` %>% str_replace("<0.01", "0.005")) %>% 
  mutate(Temp = Temperature) %>% 
  mutate(SS = `Suspended Solids` %>% str_replace("<0.5", "0.25")) %>% 
  mutate(silica = Silica %>% str_replace("<0.05", "0.025")) %>%
  mutate(Secchi_depth = `Secchi Disc Depth`) %>% 
  mutate(salinity = Salinity) %>% 
  mutate(Phaeo_pigments = `Phaeo-pigments` %>% str_replace("<0.2", "0.1")) %>%
  mutate(PH = pH) %>% 
  mutate(PO4 = `Orthophosphate Phosphorus` %>% str_replace("<0.002","0.001")) %>%
  mutate(NO2_N = `Nitrite Nitrogen` %>% str_replace("<0.002", "0.001")) %>% 
  mutate(NO3_N = `Nitrate Nitrogen` %>% str_replace("<0.002", "0.001")) %>% 
  mutate(Faecal = `Faecal Coliforms` %>% str_replace_all(c("<1" = "0.5",
                                                           "<2" = "1"))) %>%
  mutate(E_coli = `E. coli`%>% str_replace_all(c("<1" = "0.5",
                                                 "<2" = "1"))) %>%
  mutate(DO = `Dissolved Oxygen`) %>% 
  mutate(DO_sat = `Dissolved Oxygen_SAT`) %>% 
  mutate(Chl_a = `Chlorophyll-a` %>% str_replace("<0.2", "0.1")) %>%
  mutate(NH3_N = `Ammonia Nitrogen` %>% str_replace("<0.005", "0.0025")) %>%
  mutate(BOD5 = `5-day Biochemical Oxygen Demand` %>% str_replace_all(c("<0.1" = "0.05",
                                                                        ">6.5" = "13",
                                                                        "<2" = "1",
                                                                        ">8.4" = "16.8",
                                                                        ">6.3" = "12.6"))) %>% 
  select(1, 26:49) 

######################################################################
######################################################################
env1 %>% 
  mutate(Station = Sample %>% str_split("_") %>% sapply('[', 1)) %>% 
  mutate(Year = Sample %>% str_split("_") %>% sapply('[', 2)) %>% 
  mutate(Month = Sample %>% str_split("_") %>% sapply('[', 3) %>% 
           str_replace_all(c("^1$"="01",   ## some use"1,2,3" while some use "01.02.03")
                             "^2$"="02",
                             "^3$"="03",
                             "^4$"="04",
                             "^5$"="05",
                             "^6$"="06",
                             "^7$"="07",
                             "^8$"="08",
                             "^9$"="09"))) %>% 
  mutate(YM = paste0(Year, Month)) %>% 
  filter(Station %in% c("SM6", "SM17", "SM19", "MM4","MM14","MM16")) %>%
  filter(YM %in% c("201810","201811","201812","201901","201902","201903","201904",
                   "201905","201906","201907","201908","201909")) %>% 
  mutate(Season = Month %>% str_replace_all(c("12" = "Winter",
                                              "01" = "Winter",
                                              "02" = "Winter",
                                              "03" = "Spring",
                                              "04" = "Spring",
                                              "05" = "Spring",
                                              "06" = "Summer",
                                              "07" = "Summer",
                                              "08" = "Summer",
                                              "09" = "Fall",
                                              "10" = "Fall",
                                              "11" = "Fall"))) %>% 
  mutate(Region = Station %>% str_replace_all(c("SM17" = "Southern",
                                                "SM6" = "Southern",
                                                "SM19" = "Southern",
                                                "MM4" = "Eastern",
                                                "MM14" = "Eastern",
                                                "MM16" = "Eastern"))) %>% 
  select(Station,Region,Season,2:25) %>% 
  rename(Salinity = salinity) %>% 
  rename(pH = PH)-> env2

env2 %>% colnames()
df = env2[,-c(1:3)]
df = df %>% apply(2, as.numeric)
pca_res <- prcomp(df, scale = T)  

######################################################################
######################################################################
library(ggfortify)
autoplot(pca_res, data = env2, colour = 'Region', shape = 'Season',size = 5,
         loadings = TRUE, loadings.colour = 'cornflowerblue', loadings.label.colour = "Black",
         loadings.label = TRUE, loadings.label.size = 5,
         loadings.arrow = grid::arrow(length = grid::unit(8, "points")))+
  geom_hline(yintercept = 0, linetype = 3, size = 1)+
  geom_vline(xintercept = 0, linetype = 3, size = 1)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title = element_text(size = 20))+ 
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20),
        axis.title=element_text(size=20))

######################################################################
### ANOVA tests seasonal and spatial effects on environmental variables
######################################################################
# Eastern region
env_east <-
  env2 %>% 
  filter(Region == "Eastern")

df = env_east[,-c(1:3)]
df = df %>% apply(2, as.numeric)
pca_res_e <- prcomp(df, scale = T)  

pca_res_e$x[,1]

env_east %>% 
  select(Station, Season) %>% 
  mutate(PC1 = pca_res_e$x[,1]) -> dff

anova_e = aov(PC1~Station+ Season, data = dff)
summary(anova_e) # Station: F=0.04, p=0.99; Season: F=4.007, p=0.016

# Southern region
env_s <-
  env2 %>% 
  filter(Region == "Southern")

df = env_s[,-c(1:3)]
df = df %>% apply(2, as.numeric)
pca_res_s <- prcomp(df, scale = T)  

pca_res_s$x[,1]

env_s %>% 
  select(Station, Season) %>% 
  mutate(PC1 = pca_res_s$x[,1]) -> dff

anova_s = aov(PC1~Station+ Season, data = dff)
summary(anova_s) # Station: F=2.54, p=0.1; Season: F=18.7, p<0.001
