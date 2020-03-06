# Packages ----

library("addinslist")
library(tidyverse)
library(readr)
library(RColorBrewer)
library(modelr)
library(kableExtra)
library(ggThemeAssist)
library(stringr)


getwd()

# LS = Lipid Search Software i.e. table of lipids etc obtained from LipidSearch
# CD =  Compund Discoverer

##%######################################################%##
#                                                          #
####                   1st function                     ####
#                                                          #
##%######################################################%##


#  load and transorm the data obtained in LipidSearch in a format can work on in R 

readLipidSearch <- function(filename) {
  skiprow <- read_tsv(filename)
  skiprow <- which(grepl("Result", skiprow$`Job Name`)) +3
  df <- read_tsv(filename, skip = skiprow) %>%
    separate(`FA Group Key`, c("TC", "DB"), sep = ":", remove = FALSE)
  }


LSdata <- readLipidSearch("./Data/LS_TableDomBatch01.txt")
#saved in Z:\Documents\Lipidomics Domestication Batch01\R\Data


ggplot(LSdata, aes(BaseRt, `Calc Mass`, colour = Class)) +geom_point() 

# first check to see if some compounds look a bit "dodgy"; so separate the different lipid classes
# but want to look to see if there are any difference

ggplot(filter(LSdata, Class == "TG"), aes(BaseRt, `Calc Mass`, colour = DB)) +
  geom_point() +
  facet_wrap(~Class)


##%######################################################%##
#                                                          #
####                    2nd function                    ####
#                                                          #
##%######################################################%##
 

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



##%######################################################%##
#                                                          #
####                      3rd part                      ####
#                                                          #
##%######################################################%##


# now that CD has run we want to plot this data and see what makes sense
#    !!! replace Batch01 with Batch02 accordingly in the following script !!!!!!!

CDdata <- read_tsv("Data/CompDiscoBatch02Results12Dec.csv", col_types = cols(Name = col_character()))

key <- read_csv("Data/Table_CD11Dec.csv", col_types = cols(Name = col_character())) %>%  
  select(Name, Class, FA, `FA Group Key`) %>% 
  unique()
  # rename(LSFormula= Formula) %>% 
  # mutate(`FA Group Key` = str_remove(str_replace(`FA Group Key`, ":", "_"), ":00"))
   

CDdata_gathered06B2 <- CDdata %>% 
  mutate(Name = gsub("Oleic acid", "FA(18:1)", x = Name)) %>% ###use quotes not the other thing (back thingy``)
  gather(Sample, Area, contains("Area: ")) %>%
  mutate(Sample = str_replace(Sample, regex("\\-", ignore_case = TRUE), "_")) %>% 
  group_by(Name, Formula, `Molecular Weight`, `RT [min]`) %>%
  filter(!is.na(Area),
         `RT [min]` < 25,
         !is.na(Name)) %>% 
  mutate(count = n()) %>% 
  left_join(key) %>% 
  mutate(Bond = str_extract(`FA Group Key`, "[a-zA-Z]+")) %>% 
  mutate(TC = word(`FA Group Key`, 1, 1, sep = ":"),
         DB = word(`FA Group Key`, -1, -1, sep = ":")) %>% 
 # mutate(DB = as.numeric(str_extract(DB, "\\d+"))) %>%   #this separate the letters of plasmalogen from the DB IF run
  group_by(Sample) %>%
  mutate(Area2 = Area/sum(Area) * 1000000) %>% 
  ungroup() %>% 
  group_by(Formula, `Molecular Weight`, `RT [min]`) %>%
  mutate(reference = median(Area2),
         quotient = Area2/reference) %>%
  ungroup() %>%
  group_by(Sample) %>%
  mutate(quotient.median = median(quotient),
         pqn = Area2/quotient.median) # %>%
  # mutate (pop = (str_extract(Sample, "\\w+_\\w+"))) %>% 
  # separate(pop, c("population", "sample"), sep="(?<=[A-Za-z])(?=[0-9])")
  
CDdata_gathered06B2 
write_csv(CDdata_gathered06B2, path ="Data/CDdata_gathered06B2.csv")


#lipids.to.keep <- c("DG", "TG", "PC", "PI", "PA", "PG", "PS", "SM", "LPC", "FA","PE", "Cer")

res <-  CDdata_gathered06B2 %>% 
  #filter(grepl(pattern = paste(lipids.to.keep, collapse = "|"), Name)) %>% 
  filter(!grepl("Similar", Name))
           
Batch02 <- res %>% 
  select(Name, Formula, Class, `FA Group Key`, TC, DB, count, `Molecular Weight`, Sample, `RT [min]`, Area, pqn, Area2)

write_csv(Batch02, path = "Data/Batch02.csv")
view(Batch02)


# Pre-final adjustment (still need to remove NA)

p1B2 <- Batch02 %>% 
  group_by(Class) %>% 
  rename( `Lipid species` = `FA Group Key`)  
 # filter(!grepl("e|p|t|Q|d|O", `Lipid species`)) #This removes letters from the Lipid Species and thus DB IF run

write_csv(p1B2, path = "Data/p1B2.csv")

# PLOTS: 
# First plot with areas, not conc yet;
# this gives the boxplot w title "Batch02", not a nice plot, only good for a 1st check

batch02.plot <-  ggplot(p1B2, aes(`Lipid species`, log(Area))) +
  geom_boxplot(aes(fill= `Lipid species`)) 
  #facet_wrap(~Class, scales ="free")

batch02.plot +
  theme(legend.position ="none", 
        axis.text.x = element_text(angle = 45, vjust =0.5)) +
  labs( title = paste("Batch02"))

#### plot normalised data 

batch02.plot.pqn <-  ggplot(p1B2, aes(`Class`, log(Area2))) +
  geom_boxplot(aes(fill= `Lipid species`)) 
  #facet_wrap(~Class, scales ="free")

batch02.plot.pqn +
  theme(#legend.position ="none", 
        axis.text.x = element_text(angle = 45, vjust =0.5)) +
  labs( title = paste("Batch02 - normalised data"))

# Need to add populations  by Joining the table (Batch02Standards) containing
# the STANDARDS for Batch01   (run with the compounds in the LCMS)                
#                 with  Batch02Results i.e. to the compunds of Batch02                 
#                                                          #
##%######################################################%##

Batch02Standards <- read_csv("./Data/Batch02Standards.csv")
#view(Batch02Standards)

batch02.standards <- Batch02Standards %>% 
  select(Class, StandardsName, Conc, slope, intercept)

# join the 2 tables, in order to add the populations 

batch02.values <-  p1B2 %>% 
  left_join(batch02.standards, by = "Class")  %>% 
  mutate(conc_mgmL_compounds = as.numeric(exp((log(Area) - intercept)/ slope))) %>% 
  select(-slope, - intercept) %>% 
  #filter(!Class %in% c("LPE", "LPI", "LPS", "PA")) %>% 
  filter(!is.na(Class)) 

#  this below excludes the plasmalogen but not the Ceramides and the SPHM (kept by deleting d and t from grepl)

batch02.values.noplasm <-  batch02.values  %>% 
    filter(!grepl("e|p|Q|O", `Lipid species`))

# split the Sample column in population w domestication and sample number

batch02.pop <- batch02.values.noplasm %>% 
        mutate (pop = (str_extract(Sample, "\\w+_\\w+"))) %>% 
        separate(pop, c("population", "sample"), sep="(?<=[A-Za-z])(?=[0-9])")
  
 
# this plot (still ugly) is based on areas not on concentrations (although it is in the name, please disregard it)

batch01.plot.conc <-  ggplot(batch02.pop, aes(`Class`, log(pqn))) + 
  geom2boxplot(aes(fill= `population`), position = position_dodge(width=0.75)) +
  #facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1), 
    axis.text = element_text(size = 7, vjust = -0.25), 
    axis.text.x = element_text(angle = 40), 
   # panel.background = element_rect(fill = "palegoldenrod"), 
    legend.position = "right") +
    labs(title = "Batch02: lipid classes", 
    y = "log2 norm. area") 

batch02.plot.conc


#write_csv(batch01.plot.conc.pqn, path = "Data/batch01.plot.conc.pqn.csv")

# sum the areas belonging to the same class of lipids (i.e. the areas of the different species)

batch02.sum <-  batch02.pop %>% 
  group_by(Class, sample) %>%  
  mutate(sumpqn = sum(pqn)) %>% 
  mutate(population = tolower(population))
         
         
 unique(batch02.sum$population)        

  
#sum_conc_mgml = sum(conc_mgmL_compounds))

batch02.sum #this is THE data for statistical analysis of ALL lipids based on normalised areas
write_csv(batch02.sum, path = "Data/batch02.sum.csv") # file saved :)

# NICE plot finally!

batch02.plot.all <-  ggplot(batch02.sum, aes(`Class`, log(sumpqn)))+
  geom_boxplot(aes(fill = `population`), width =1) +
  scale_fill_brewer(palette = "Paired") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch 02 Lipid Classes", y = "Log normalised area", x = "Lipid Classes") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11)) 
        #theme(legend.position = c(0.12, 0.9), legend.direction = "horizontal") 

batch02.plot.all







