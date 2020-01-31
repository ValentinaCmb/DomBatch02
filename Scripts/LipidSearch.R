library(tidyverse)
library(readr)
getwd()

# LS = Lipid Search Software i.e. table of lipids etc obtained from LipidSearch
# CD =  Compund Discoverer

## ---- 1st function -----------

# to transorm the data obtained in lipid search in a format that is usable/we can work on

readLipidSearch <- function(filename) {
  skiprow <- read_tsv(filename)
  skiprow <- which(grepl("Result", skiprow$`Job Name`)) +3
  df <- read_tsv(filename, skip = skiprow) %>%
    separate(`FA Group Key`, c("TC", "DB"), sep = ":", remove = FALSE)
  }


LSdata <- readLipidSearch("./Data/LS_TableDomBatch01.txt")


ggplot(LSdata, aes(BaseRt, `Calc Mass`, colour = Class)) +geom_point() 

#some compounds look a bit dodgy so lets separate the different lipid classes
#but want to look to see if there are any difference
ggplot(filter(LSdata, Class == "TG"), aes(BaseRt, `Calc Mass`, colour = DB)) +geom_point() +
  facet_wrap(~Class)


# ---------- 2nd function ------------------------------------------------------------------

# following step 1, now we want to pull out  data from the lipidsearch table to make a new table that will run in CD

LS_to_CD_list <- function(df){
  #unassigned <- filter(df, Rej == 0) #need to change to do this only if there is a column called Rej
  calc <- df %>%
    ungroup() %>%
    select(Name = LipidMolec, LSFormula = Formula, RT = BaseRt, MainIon, `Calc Mass`, Class, FA, `FA Group Key`) %>%
    unique() %>%
    mutate(C = as.numeric(str_match(LSFormula, "\\bC(\\d+)")[,2]),
           H = as.numeric(str_match(LSFormula, "\\bH(\\d+)")[,2]),
           N = as.numeric(ifelse(grepl("N", LSFormula), 
                                 ifelse(grepl("N\\d+", LSFormula), 
                                        str_match(LSFormula, "\\bN(\\d+)")[,2], 1), 0)), 
           O = as.numeric(ifelse(grepl("O", LSFormula), 
                                 ifelse(grepl("O\\d+", LSFormula), 
                                        str_match(LSFormula, "\\bO(\\d+)")[,2], 1), 0)),
           P = as.numeric(ifelse(grepl("P", LSFormula), 
                                 ifelse(grepl("P\\d+", LSFormula), 
                                        str_match(LSFormula, "\\bP(\\d+)")[,2], 1), 0)),
           S = ifelse(grepl("S", LSFormula), 
                      ifelse(grepl("S\\d+", LSFormula), 
                             str_match(LSFormula, "\\bS(\\d+)")[,2], 1), ""),
           N = ifelse(MainIon == "+NH4", N + 1, N),
           H = ifelse(MainIon == "+NH4", H + 3, H),
           Name = ifelse(MainIon == "+NH4", paste(Name, "+NH4", sep = " "), Name),
           Formula = paste("C", C, "H", H, ifelse(N > 0, "N", ""), ifelse(N > 1, N, ""), 
                           "O", O, ifelse(P > 0, "P", ""), ifelse(P > 1, P, ""), sep = "")) 
  # %>%
    # mutate(Lipid_Class = ifelse(grepl("NO6", Formula) & !grepl("P", Formula), "TG", 
                               # ifelse(grepl("N2O5P", Formula), "CAEP", NA)),
           # Lipid_Class = ifelse(grepl("NO8P", Formula) & C %% 2 == 0, "PC" ,
               #                 ifelse(grepl("NO8P", Formula) & C %% 2 != 0, "PE", Lipid_Class)),
           # Lipid_Class = ifelse(C>30 & grepl("NO7P", Formula) & C %% 2 == 0, "PC",
                              #  ifelse(grepl("NO7P", Formula) & C %% 2 != 0, "PE", Lipid_Class)),
          # Lipid_Species = ifelse(grepl("NO6", Formula) & !grepl("P", Formula),
                                 # paste(C-3, ((C-3)*2-(H-5))/2, sep = ":"), NA)) 
  
  
  return(calc)
}

LSdata01 <- LS_to_CD_list(LSdata)
write_csv(LSdata01, path = "Data/Table_CD28Nov.csv")

# this below is not needed 
T01 <- LSdata01 %>% 
      select(Name, LSFormula, RT, MainIon, `Calc Mass`, Class, FA, `FA Group Key`) %>% 
      filter(Class %in% c("DG", "FA", "PE", "LPC", "LPE", "LPI", "LPS", "PA", "PC", "PG", "PI", "PS", "TG"))



# ------------------------ 3rd part ------------------------------------------------- #

# now that CD has run we want to plot this data and see what makes sense


require(tidyverse)
require(RColorBrewer)

CDdata <- read_tsv("Data/CompDiscoBatch01Results11Dec.csv", col_types = cols(Name = col_character()))

key <- read_csv("Data/Table_CD28Nov.csv", col_types = cols(Name = col_character())) %>%    
  rename(LSFormula= Formula) %>% 
  mutate(`FA Group Key` = str_remove(str_replace(`FA Group Key`, ":", "_"), ":00"))
  

CDdata_gathered06 <- CDdata %>% 
  gather(Sample, Area, contains("Area: ")) %>%
  group_by(Name, Formula, `Molecular Weight`, `RT [min]`) %>%
  filter(!is.na(Area),
         `RT [min]` < 25,
         !is.na(Name)) %>% 
  mutate(count = n()) %>% 
  left_join(key) %>% 
  mutate(TC = word(`FA Group Key`, 1, 1, sep = "_"),
         DB = word(`FA Group Key`, -1, -1, sep = "_"))  
#mutate(Batch = str_extract(sample, "[B]\\d+"),
#       StandardsName =str_extract(sample, "\\d+"))

lipids.to.keep <- c("DG", "TG", "PC", "PI", "PA", "PG", "PS", "SM", "LPC", "Oleic acid","PE", "Cer")

res <-  CDdata_gathered06 %>% 
  filter(grepl(pattern = paste(lipids.to.keep, collapse = "|"), Name)) %>% 
  filter(!grepl("Similar", Name))
           
Batch01 <- res %>% 
  select(Name, Formula, Class, TC, DB, count, `Molecular Weight`, Sample, `RT [min]`, Area)
  


PC <- CDdata_gathered %>% 
  filter(Class == "PC")

ggplot( aes(`RT [min]`, `Molecular Weight`, colour = DB)) +
  geom_point()

p1 <- CDdata_gathered %>% 
  group_by(Batch)

ggplot(p1, aes(`StandardsName`, log(area), colour = Batch)) +
  #geom_boxplot()+
  geom_point() +
  facet_wrap(~Class, scales ="free")


# ----------------previous function------------------------------------------------------------------ #

CDdata_gathered03 <- CDdata %>% 
  gather(Sample, Area, contains("Area: ")) %>%
  group_by(Name, Formula, `Molecular Weight`, `RT [min]`) %>%
  filter(!is.na(Area),
         `RT [min]` < 25
  ) %>% 
  mutate(count = n()) %>% 
  left_join(key) %>% 
  mutate(TC = word(`FA Group Key`, 1, 1, sep = "_"),
         DB = word(`FA Group Key`, -1, -1, sep = "_"))  
#mutate(Batch = str_extract(sample, "[B]\\d+"),
#       StandardsName =str_extract(sample, "\\d+"))

