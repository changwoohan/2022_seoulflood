rm(list=ls())

# load the needed packages and functions
project.folder = paste0(print(here::here()),'/')
source(paste0(project.folder,'/code/packages/packages_to_load_20231205HCW.R'))
source(paste0(project.folder,'/code/functions/functions_20231227.R'))

#########################################################################################################
# data load: main analysis 
dir(here::here("data/raw_data"))

name_1  <- paste0("data/raw_data/seflood_20231226.csv")
prev    <- read.csv(here::here(name_1), fileEncoding = 'cp949') 


name_1  <- paste0("data/other_data/disease_name_20231226HCW.xlsx")
d_name  <- read_excel(here::here(name_1), sheet=1)
d_name$z <- gsub(";", "\n", d_name$z )


# define the intervention

sm_n <- prev %>%  
  filter(hit %in% c(1,2,3)) %>% 
  mutate(eq =ifelse(DT<as.Date("2022-08-08"), 0, 1),
         intt=ifelse(hit==1,0,1),
         int =ifelse(intt==1 & eq==1, 1, 0)) 

s_n <- prev %>%  
  filter(hit %in% c(1,3)) %>% 
  mutate(eq =ifelse(DT<as.Date("2022-08-08"), 0, 1),
         intt=ifelse(hit==3,1,0),
         int =ifelse(intt==1 & eq==1, 1, 0)) 
  
m_n <- prev %>%  
  filter(hit %in% c(1,2)) %>% 
  mutate(eq =ifelse(DT<as.Date("2022-08-08"), 0, 1),
         intt=ifelse(hit==2,1,0),
         int =ifelse(intt==1 & eq==1, 1, 0)) 


# period from -16 to +11
all_code <- colnames(prev)[str_detect(colnames(prev), c("ALL_")) 
                     & str_detect(colnames(prev), c("_T"))]

spe_code <- c("A_1_T", "A_4_T", "B_4_T",  "B_5_T",
              "F_1_T",  "F_2_T",   "F_3_T",   "F_4_T",   "F_5_T",
              "J_1_T",   "J_2_T",   "J_3_T",  "J_4_T",   "J_5_T",
              "E_1_T",   "E_2_T",   "E_3_T",   "E_4_T",
              "S_1_T", "S_2_T",   "S_3_T",   "S_4_T",   "S_5_T",
              "S_6_T",   "S_7_T",   "S_8_T",   "S_9_T",  "S_10_T")

# run the gsynth analysis 
han_gsynth_95(dat.use=sm_n, val.name=all_code, date="main", fol="ALL_SM_", fig_name="sm_n")
han_gsynth_95(dat.use=s_n,  val.name=all_code, date="main", fol="ALL_SEVERE_", fig_name="s_n")
han_gsynth_95(dat.use=m_n,  val.name=all_code, date="main", fol="ALL_MILD_",   fig_name="m_n")

han_gsynth_95(dat.use=sm_n, val.name=spe_code, date="main", fol="SPE_SM_", fig_name="spe_sm_n")
han_gsynth_95(dat.use=s_n,  val.name=spe_code, date="main", fol="SPE_SEVERE_", fig_name="spe_s_n")
han_gsynth_95(dat.use=m_n,  val.name=spe_code, date="main", fol="SPE_MILD_",   fig_name="spe_m_n")


han_gsynth_ALL_fig(dat.use=sm_n,  folder="fin_main", fig_name="FigureS2.Main_ALL_SM")
han_gsynth_ALL_fig(dat.use=s_n,   folder="fin_main", fig_name="Figure2.Main_ALL_severe")
han_gsynth_ALL_fig(dat.use=m_n,   folder="fin_main", fig_name="FigureS3.Main_ALL_mild")


han_gsynth_SPE_fig(dat.use=s_n,  folder="fin_main", fig_name="Figure3.Main_SPE_severe")
han_gsynth_SPE_fig(dat.use=m_n,  folder="fin_main", fig_name="FigureS4.Main_SPE_mild")



#########################################################################################################
# data load: sensitivity 1
dir(here::here("data/raw_data"))

name_1  <- paste0("data/raw_data/seflood_20231226_sensi_20.csv")
prev    <- read.csv(here::here(name_1), fileEncoding = 'cp949') 

name_1  <- paste0("data/other_data/disease_name_20231226HCW.xlsx")
d_name  <- read_excel(here::here(name_1), sheet=1)
d_name$z <- gsub(";", "\n", d_name$z )


# define the intervention

s_n_s1 <- prev %>%  
  filter(hit %in% c(1,3)) %>% 
  mutate(eq =ifelse(DT<as.Date("2022-08-08"), 0, 1),
         intt=ifelse(hit==3,1,0),
         int =ifelse(intt==1 & eq==1, 1, 0)) 

m_n_s1 <- prev %>%  
  filter(hit %in% c(1,2)) %>% 
  mutate(eq =ifelse(DT<as.Date("2022-08-08"), 0, 1),
         intt=ifelse(hit==2,1,0),
         int =ifelse(intt==1 & eq==1, 1, 0)) 


all_code <- colnames(prev)[str_detect(colnames(prev), c("ALL_")) 
                           & str_detect(colnames(prev), c("_T"))]

spe_code <- c("A_1_T", "A_4_T", "B_4_T",  "B_5_T",
              "F_1_T",  "F_2_T",   "F_3_T",   "F_4_T",   "F_5_T",
              "J_1_T",   "J_2_T",   "J_3_T",  "J_4_T",   "J_5_T",
              "E_1_T",   "E_2_T",   "E_3_T",   "E_4_T",
              "S_1_T", "S_2_T",   "S_3_T",   "S_4_T",   "S_5_T",
              "S_6_T",   "S_7_T",   "S_8_T",   "S_9_T",  "S_10_T")

# run the gsynth analysis 
han_gsynth_95(dat.use=s_n_s1,  val.name=all_code, date="sens1", fol="ALL_SEVERE_", fig_name="s_n")
han_gsynth_95(dat.use=m_n_s1,  val.name=all_code, date="sens1", fol="ALL_MILD_",   fig_name="m_n")

han_gsynth_95(dat.use=s_n_s1,  val.name=spe_code, date="sens1", fol="SPE_SEVERE_", fig_name="s_n")
han_gsynth_95(dat.use=m_n_s1,  val.name=spe_code, date="sens1", fol="SPE_MILD_", fig_name="m_n")




#########################################################################################################
# data load: sensitivity 2
dir(here::here("data/raw_data"))

name_1  <- paste0("data/raw_data/seflood_20231226_sensi_falsi.csv")
prev    <- read.csv(here::here(name_1), fileEncoding = 'cp949') 

name_1  <- paste0("data/other_data/disease_name_20231226HCW.xlsx")
d_name  <- read_excel(here::here(name_1), sheet=1)
d_name$z <- gsub(";", "\n", d_name$z )


# define the intervention

s_n_s2 <- prev %>%  
  filter(hit %in% c(1,3)) %>% 
  mutate(eq =ifelse(DT<as.Date("2021-08-08"), 0, 1),
         intt=ifelse(hit==3,1,0),
         int =ifelse(intt==1 & eq==1, 1, 0)) 

m_n_s2 <- prev %>%  
  filter(hit %in% c(1,2)) %>% 
  mutate(eq =ifelse(DT<as.Date("2021-08-08"), 0, 1),
         intt=ifelse(hit==2,1,0),
         int =ifelse(intt==1 & eq==1, 1, 0)) 


all_code <- colnames(prev)[str_detect(colnames(prev), c("ALL_")) 
                           & str_detect(colnames(prev), c("_T"))]

spe_code <- c("A_1_T", "A_4_T", "B_4_T",  "B_5_T",
              "F_1_T",  "F_2_T",   "F_3_T",   "F_4_T",   "F_5_T",
              "J_1_T",   "J_2_T",   "J_3_T",  "J_4_T",   "J_5_T",
              "E_1_T",   "E_2_T",   "E_3_T",   "E_4_T",
              "S_1_T", "S_2_T",   "S_3_T",   "S_4_T",   "S_5_T",
              "S_6_T",   "S_7_T",   "S_8_T",   "S_9_T",  "S_10_T")

# run the gsynth analysis 
han_gsynth_95(dat.use=s_n_s2,  val.name=all_code, date="sens2", fol="ALL_SEVERE_", fig_name="s_n")
han_gsynth_95(dat.use=m_n_s2,  val.name=all_code, date="sens2", fol="ALL_MILD_",   fig_name="m_n")

han_gsynth_95(dat.use=s_n_s2,  val.name=spe_code, date="sens2", fol="SPE_SEVERE_", fig_name="s_n")
han_gsynth_95(dat.use=m_n_s2,  val.name=spe_code, date="sens2", fol="SPE_MILD_", fig_name="m_n")










