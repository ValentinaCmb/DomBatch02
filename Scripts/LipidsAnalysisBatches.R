# Packages ----

#Load require packages

library("addinslist")
library(tidyverse)
library(readr)
library(RColorBrewer)
library(modelr)
library(kableExtra)
library(ggThemeAssist)



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


CDdata <- read_tsv("Data/CompDiscoBatch02Results12Dec.csv", col_types = cols(Name = col_character()))

key <- read_csv("Data/Table_CD11Dec.csv", col_types = cols(Name = col_character())) %>%  
  select(Name, Class, FA, `FA Group Key`) %>% 
  unique()
  # rename(LSFormula= Formula) %>% 
  # mutate(`FA Group Key` = str_remove(str_replace(`FA Group Key`, ":", "_"), ":00"))
  # 

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
  mutate(DB = as.numeric(str_extract(DB, "\\d+"))) %>% 
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
  # 
 
CDdata_gathered06B2
write_csv(CDdata_gathered06B2, path ="Data/CDdata_gathered06B2.csv")


lipids.to.keep <- c("DG", "TG", "PC", "PI", "PA", "PG", "PS", "SM", "LPC", "FA","PE", "Cer")

res <-  CDdata_gathered06B2 %>% 
  filter(grepl(pattern = paste(lipids.to.keep, collapse = "|"), Name)) %>% 
  filter(!grepl("Similar", Name))
           
Batch02.selected <- res %>% 
  select(Name, Formula, Class, `FA Group Key`, TC, DB, count, `Molecular Weight`, Sample, `RT [min]`, Area, pqn, Area2)

write_csv(Batch02.selected, path = "Data/Batch02.selected.csv")
view(Batch02.selected)
unique(Batch02.selected$Class)
# Pre-final adjustment (still need to remove NA)

p1.selected <- Batch02.selected %>% 
  group_by(Class) %>% 
  rename( `Lipid species` = `FA Group Key`) %>% 
  filter(!grepl("e|p|t|Q|d|O", `Lipid species`))

write_csv(p1.selected, path = "Data/p1.selected.csv")

# PLOTS: 
# First plot with areas, not conc yet ##
# this gives the boxplot w title "Batch01" 

batch02selected.plot <-  ggplot(p1.selected, aes(`Lipid species`, log(Area))) +
  geom_boxplot(aes(fill= `Lipid species`)) +
  facet_wrap(~Class, scales ="free")

batch02selected.plot +
  theme(legend.position ="none", 
        axis.text.x = element_text(angle = 45, vjust =0.5)) +
  labs( title = paste("Batch02"))

#### plot normalised data 
batch02selected.plot.pqn <-  ggplot(p1.selected, aes(`Lipid species`, log(pqn))) +
  geom_boxplot(aes(fill= `Lipid species`)) +
  facet_wrap(~Class, scales ="free")

batch02selected.plot.pqn +
  theme(legend.position ="none", 
        axis.text.x = element_text(angle = 45, vjust =0.5)) +
  labs( title = paste("Batch02 - normalised data"))
                                                          
##%######################################################%##
#                                                          #
#    Join the table (Batch01StandardsNew) containing      #
#       the STANDARDS for Batch01                          #
#        (ran with the compounds in the LCMS)              #  
#                 with  Batch01Results                     #
#          i.e. to the compunds of Batch01                 #
#                                                          #
##%######################################################%##

Batch02Standards <- read_csv("./Data/Batch02Standards.csv")
view(Batch02Standards)

batch02.standards.selected <- Batch02Standards %>% 
  select(Class, StandardsName, Conc, slope, intercept)


# calculate the concentrations of the compounds in Batch01 based on the standards in Batch01

batch02.values.selected <-  p1.selected %>% 
  left_join(batch02.standards.selected, by = "Class")  %>% 
  mutate(conc_mgmL_compounds = as.numeric(exp((log(Area) - intercept)/ slope))) %>% 
  select(-slope, - intercept) %>% 
  filter(!Class %in% c("LPE", "LPI", "LPS", "PA")) %>% 
  filter(!is.na(Class))

#check elimination of those lipids with unique
unique(batch02.values.selected$Class)

batch02.values.selected


# split the Sample column in population w domestication and sample number

batch02.selected.pop <- batch02.values.selected %>% 
        mutate (pop = (str_extract(Sample, "\\w+_\\w+"))) %>% 
        separate(pop, c("population", "sample"), sep="(?<=[A-Za-z])(?=[0-9])") %>% 
  mutate(population = tolower(population))
  
 unique(batch02.selected.pop$Class)

#str_extract("AB:C_N01.DE", "\\w+_\\w+" )
#unlist(str_split("pop","(?<=[A-Za-z])(?=[0-9])"))
# str_extract("AB:C_N01.DE", "\\d+" )
# str_extract("C_N01", "\\d+" )
# str_remove_all("C_N01", "[[:punct:]]" )
# str_split(str_split("C_N01","_")[[1]][[2]], "(?<=[A-Za-z])(?=[0-9])")
# unlist(str_split("C_N01","(?<=[A-Za-z])(?=[0-9])"))[2]


##%######################################################%##
#                                                          #
#               PLOTS of the CONCENTRATIONS             ####
#                                                         #
##%######################################################%##

# ! batch02.pop has the population now

batch02.selected.plot.conc <-  ggplot(batch02.selected.pop, aes(`Lipid species`, log(conc_mgmL_compounds))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1), 
    axis.text = element_text(size = 7, vjust = -0.25), 
    axis.text.x = element_text(angle = 40), 
    panel.background = element_rect(fill = "palegoldenrod"), 
    legend.position = "right") +labs(title = "Batch02: lipid species concentrations", 
    y = "log concentrations") 

batch02.selected.plot.conc

###normalised area results -  NOT for concentration

batch02.selected.plot.conc.pqn <-  ggplot(batch02.selected.pop, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch02: lipid species normalised areas", 
                                         y = "log(pqn)") 

batch02.selected.plot.conc.pqn
write_csv(batch02.selected.pop, path = "Data/batch02.selected.pop.csv")

# SINGLE LIPIDS CLASSES

# NO FA present
#FAnorm.area <- batch02.selected.pop %>% 
#  filter(Class == "FA")

#FA.batch02.plot.conc.pqn <-  ggplot(FAnorm.area, aes(`Lipid species`, log(pqn))) + 
#  geom_boxplot(aes(fill= `population`)) +
  #facet_wrap(~Class, scales ="free") + 
#  theme(plot.subtitle = element_text(vjust = 1), 
#        plot.caption = element_text(vjust = 1), 
#        axis.text = element_text(size = 7, vjust = -0.25), 
#        axis.text.x = element_text(angle = 40), 
        #panel.background = element_rect(fill = "palegoldenrod"), 
#        legend.position = "right") +labs(title = "Batch02: FA lipid species normalised areas", 
#                                         y = "log(pqn)") 

#FA.batch02.plot.conc.pqn


PGnorm.area <- batch02.selected.pop %>% 
  filter(Class == "PG")

PG.batch02.plot.conc.pqn <-  ggplot(PGnorm.area, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        #panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch02: PG lipid species normalised areas", 
                                         y = "log(pqn)") 

PG.batch02.plot.conc.pqn


TGnorm.area <- batch02.selected.pop %>% 
  filter(Class == "TG")

TG.batch02.plot.conc.pqn <-  ggplot(TGnorm.area, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        #panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch02: TG lipid species normalised areas", 
                                         y = "log(pqn)") 

TG.batch02.plot.conc.pqn


DGnorm.area <- batch02.selected.pop %>% 
  filter(Class == "DG")

DG.batch02.plot.conc.pqn <-  ggplot(DGnorm.area, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        #panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch02: DG lipid species normalised areas", 
                                         y = "log(pqn)") 

DG.batch02.plot.conc.pqn


LPCnorm.area <- batch02.selected.pop %>% 
  filter(Class == "LPC")

LPC.batch02.plot.conc.pqn <-  ggplot(LPCnorm.area, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        #panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch02: LPC lipid species normalised areas", 
                                         y = "log(pqn)") 

LPC.batch02.plot.conc.pqn


PCnorm.area <- batch02.selected.pop %>% 
  filter(Class == "PC")

PC.batch02.plot.conc.pqn <-  ggplot(PCnorm.area, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        #panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch02: PC lipid species normalised areas", 
                                         y = "log(pqn)") 

PC.batch02.plot.conc.pqn


PEnorm.area <- batch02.selected.pop %>% 
  filter(Class == "PE")

PE.batch02.plot.conc.pqn <-  ggplot(PEnorm.area, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        #panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch02: PE lipid species normalised areas", 
                                         y = "log(pqn)") 

PE.batch02.plot.conc.pqn


PInorm.area <- batch02.selected.pop %>% 
  filter(Class == "PI")

PI.batch02.plot.conc.pqn <-  ggplot(PInorm.area, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        #panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch02: PI lipid species normalised areas", 
                                         y = "log(pqn)") 

PI.batch02.plot.conc.pqn


PSnorm.area <- batch02.selected.pop %>% 
  filter(Class == "PS")

PS.batch02.plot.conc.pqn <-  ggplot(PSnorm.area, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        #panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch02: PS lipid species normalised areas", 
                                         y = "log(pqn)") 

PS.batch02.plot.conc.pqn

#######
# Save the plots


ggsave("batch01.plot.FA.png", width =10, height = 8)


ggsave("batch01.plot.DG.png", width =10, height = 8)




