# Refresher on R

# install.packages("tidyverse")
# load required library
library(tidyverse)
library(haven)

# import data
mortality <- read_dta("data/mortality17.dta")
mortality <- as_factor(mortality)
mortality
mortality$sex <- as_factor(mortality$sex )
# inspect the data
View(mortality)
glimpse(mortality)
library(skimr)
skim(mortality)
summary(mortality)

# Export data

# to csv
write_csv(mortality, "data/mortality.csv")

# to spss
write_sav(mortality, "data/mortality.sav")


# Data manipulation
mortality <- read_dta("data/mortality17.dta") |> as_factor()

# select for columns
select(mortality, age, sex,ethnic)
mortality |> select( age, sex,ethnic)

mortality |> select(2:5)

mortality |> select(district, age, religion)
mortality |> select(c(1,3,5))

mortality |> select(!vcode)

# to extract numeric columns
mortality |> select(where(is.numeric))

# to extract categorical columns
mortality |> select(where(is.factor))

#
mortality |> select(bmi, weight, height, everything())

# help
help("select")
?select

?filter

# summarise and groupby
summary_stat <- 
mortality |>  filter(!is.na(weight)) |> group_by(sex, agegrp) |>  
  summarize(min_weight = min(weight), 
            max_weight = max(weight),
            mean_weight = mean(weight),
            sd_wt = sd(weight),
            n = n())

write_csv(summary_stat, "data/summary_stat.csv")

# 
library(writexl)
write_xlsx(summary_stat, "data/summary_stat.xlsx")

library(flextable)

summary_stat2 <- 
summary_stat |> flextable() |> autofit()

save_as_docx(summary_stat2, path = "data/summary_stat2.docx")

# 
summary(mortality)
mortality |> count(sex,agegrp)


mortality |> 
  mutate(bmi_status = case_when(
    bmi < 18.5 ~ "Underweight",
    bmi >= 18.5 & bmi < 25 ~ "Normal weight",
    bmi >= 25 ~ "Overweight")) |> 
  select(bmi, bmi_status)


summary(mortality$age)


mortality |> 
  mutate(ht_group = ifelse(height < 170, "short", "long")) |> 
  select(height, ht_group) |> view()

# data visualization
library(ggplot2)

# scatter plot
ggplot(mortality, aes(x = age, y = diastolic)) +
  geom_point(alpha=0.7) +
  geom_smooth(method = "lm") +
  theme_bw()

# slrm
fit <- lm(mortality$diastolic ~ mortality$age)
summary(fit)

# bar chart
ggplot(mortality, aes(agegrp, fill = agegrp)) +
  geom_bar(width = 0.7) +
  guides(fill = "none") +
  theme_classic()


