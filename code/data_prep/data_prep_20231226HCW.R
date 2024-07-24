
rm(list=ls())

# load the needed packages and functions
project.folder = paste0(print(here::here()),'/')
source(paste0(project.folder,'/code/packages/packages_to_load_20231205HCW.R'))

#########################################################################################################
# prevalence data loading 
# two data sources - prevalence and incidence 
dir(here::here("data/raw_data"))
name_1  <- paste0("data/raw_data/daily_emd_prev_seflood_final.csv")
a_1 <-read.csv(here::here(name_1)) %>% 
  mutate(REGION_EMD= as.factor(REGION_EMD)) 

head(a_1)
a_1[is.na(a_1)] <-0 # change NA to 0

a_1 %>% distinct(REGION_EMD)

# total 72 sub-regions in raw data

# load the sgg coding file 
dir(here("data/other_data"))
sgg <-read_excel(here::here("data/other_data/hangan_sgg_code_20230522HCW.xlsx"), sheet=1)

se_emd <- sgg %>% mutate(REGION_SIDO = substr(ADDR_CD,1,2),
                         REGION_SGG = substr(ADDR_CD,1,5),
                         REGION_EMD = substr(ADDR_CD,1,8),
                         KOR_name = as.factor(SGG_NM)) %>% 
  filter(REGION_SGG %in% c(11590, 11620, 11650, 11560)) %>% 
  select(REGION_SGG, REGION_SIDO, KOR_name, REGION_EMD, EMD_NM) %>% 
  distinct(REGION_EMD, REGION_SIDO, EMD_NM) 
  
# total 72 sub-regions 

a_2 <- a_1 %>%
  left_join(se_emd, by=c("REGION_EMD"))

##########################################################################################################################


# load the affected EMD region information 
name2  <- paste0("data/raw_data/fin_emd_region_20230912HCW.xlsx")

aff_buld <-read_excel(here::here(name2), sheet=1) 
aff_road <-read_excel(here::here(name2), sheet=2) 

aff_total <- aff_buld %>% 
  left_join(aff_road, by="REGION_EMD") %>% 
  filter(!substr(REGION_EMD,1,5) %in% c(11680))

head(aff_total)

dt<- aff_total %>% 
  summarise_at(vars(BULD_PERCENT, RW_PERCENT),
               list(min=min, Q1=~quantile(., probs = 0.25),
                    median=median, Q3=~quantile(., probs = 0.75),
                    max=max))

mean(aff_total$BULD_PERCENT)

# select the affected region: hard mid low
aff_total_1 <- aff_total %>% 
  mutate(hit = case_when(
    BULD_PERCENT>=15                  ~ 3,
    BULD_PERCENT<15 & BULD_PERCENT>=5 ~ 2,
    BULD_PERCENT<5                    ~ 1)) %>% 
  mutate(REGION_EMD =as.character(REGION_EMD))

table(aff_total_1$hit)


a_3 <- a_2 %>% 
  left_join(aff_total_1 , by=c("REGION_EMD")) %>% 
  group_by(REGION_EMD) %>% 
  arrange(DATE, .by_group = TRUE) %>% 
  mutate(serial =row_number()) %>% 
  mutate(time=serial-161) %>% 
  mutate(DT = DATE)  %>% 
  filter(DATE>'2022-04-30')
  
write.csv(a_3, here::here("data/raw_data/seflood_20231226.csv"), fileEncoding = 'cp949')
# dataset for main analysis 

##########################################################################################################################
# different definition

# select the affected region: hard mid low
aff_total_1 <- aff_total %>% 
  mutate(hit = case_when(
    BULD_PERCENT>=20                  ~ 3,
    BULD_PERCENT<20 & BULD_PERCENT>=5 ~ 2,
    BULD_PERCENT<5                    ~ 1)) %>% 
  mutate(REGION_EMD =as.character(REGION_EMD))

table(aff_total_1$hit)


a_3 <- a_2 %>% 
  left_join(aff_total_1 , by=c("REGION_EMD")) %>% 
  group_by(REGION_EMD) %>% 
  arrange(DATE, .by_group = TRUE) %>% 
  mutate(serial =row_number()) %>% 
  mutate(time=serial-161) %>% 
  mutate(DT = DATE) %>% 
  filter(DATE>'2022-04-30')

write.csv(a_3, here::here("data/raw_data/seflood_20231226_sensi_20.csv"), fileEncoding = 'cp949')
#_sensi is a dataset with different definition of high exposure 


##########################################################################################################################
# last sensitivity analysis for 1year before the flood

dir(here::here("data/raw_data"))
name_1  <- paste0("data/raw_data/daily_emd_prev_seflood_falsi.csv")
a_1 <-read.csv(here::here(name_1)) %>% 
  mutate(REGION_EMD= as.factor(REGION_EMD)) 

head(a_1)
a_1[is.na(a_1)] <-0 # change NA to 0

a_1 %>% distinct(REGION_EMD)

# total 72 sub-regions in raw data

# load the sgg coding file 
dir(here("data/other_data"))
sgg <-read_excel(here::here("data/other_data/hangan_sgg_code_20230522HCW.xlsx"), sheet=1)

se_emd <- sgg %>% mutate(REGION_SIDO = substr(ADDR_CD,1,2),
                         REGION_SGG = substr(ADDR_CD,1,5),
                         REGION_EMD = substr(ADDR_CD,1,8),
                         KOR_name = as.factor(SGG_NM)) %>% 
  filter(REGION_SGG %in% c(11590, 11620, 11650, 11560)) %>% 
  select(REGION_SGG, REGION_SIDO, KOR_name, REGION_EMD, EMD_NM) %>% 
  distinct(REGION_EMD, REGION_SIDO, EMD_NM) 

# total 72 sub-regions 

a_2 <- a_1 %>%
  left_join(se_emd, by=c("REGION_EMD"))


# load the affected EMD region information 
name2  <- paste0("data/raw_data/fin_emd_region_20230912HCW.xlsx")

aff_buld <-read_excel(here::here(name2), sheet=1) 
aff_road <-read_excel(here::here(name2), sheet=2) 

aff_total <- aff_buld %>% 
  left_join(aff_road, by="REGION_EMD") %>% 
  filter(!substr(REGION_EMD,1,5) %in% c(11680))

# select the affected region: hard mid low
aff_total_1 <- aff_total %>% 
  mutate(hit = case_when(
    BULD_PERCENT>=15                  ~ 3,
    BULD_PERCENT<15 & BULD_PERCENT>=5 ~ 2,
    BULD_PERCENT<5                    ~ 1)) %>% 
  mutate(REGION_EMD =as.character(REGION_EMD))

table(aff_total_1$hit)

a_3 <- a_2 %>% 
  left_join(aff_total_1 , by=c("REGION_EMD")) %>% 
  group_by(REGION_EMD) %>% 
  arrange(DATE, .by_group = TRUE) %>% 
  mutate(serial =row_number()) %>% 
  mutate(time=serial-100) %>% 
  mutate(DT = DATE) %>% 
  filter(DATE>'2021-04-30')

write.csv(a_3, here::here("data/raw_data/seflood_20231226_sensi_falsi.csv"), fileEncoding = 'cp949')
# dataset for main analysis 