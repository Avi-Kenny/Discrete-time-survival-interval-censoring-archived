###############################.
##### Data quality checks #####
###############################.

# !!!!! This file needs updating

chk(2, "DQA: START")

# Setup
# !!!!! Needs to change; first_hiv_pos_dt and last_hiv_neg_dt vars removed
dqa <- function(test) { if (test==F) { stop("DQA Error") } }
dat_grp2 <- dat_raw %>% dplyr::group_by(id) %>%
  dplyr::summarize(
    unique_T_plus = (function(vec){ length(unique(vec)) })(first_hiv_pos_dt), # !!!!!
    unique_T_minus = (function(vec){ length(unique(vec)) })(last_hiv_neg_dt) # !!!!!
  )

# Tests on raw data
dqa(identical(sort(unique(dat_raw$sex)), c("Female", "Male")))
dqa(sum(dat_raw$died==0, na.rm=T) + sum(dat_raw$died==1, na.rm=T)
    == length(unique(dat_raw$id)))
dqa(max(dat_grp2$unique_T_plus)==1)
dqa(max(dat_grp2$unique_T_minus)==1)

# Tests on processed data
dqa(length(unique(dat_prc$id))==max(dat_prc$id))
dqa(sum(dat_grp$T_plus-dat_grp$T_minus<0, na.rm=T)==0)
dqa(sum(dat_grp$case==999)==0)

# Tests comparing dat_raw vs. dat_prc
dqa(sum(dat_prc$y==1, na.rm=T)==sum(dat_raw$died==1, na.rm=T))
dqa(sum(dat_prc$y==0, na.rm=T)==
      sum(dat_raw$died==0, na.rm=T) + sum(is.na(dat_raw$died)))



# Distribution of "year at entry"
dat_prc %>%
  group_by(id) %>%
  mutate(min_year=min(year)) %>%
  filter(year==min_year) %>%
  xtabs(~year, data=.)

# Distribution of "year of test"
xtabs(~ResultDate, data=dat_prc)

# Distribution of "year of first test"
dat_prc %>%
  group_by(id) %>%
  filter(!is.na(ResultDate)) %>%
  mutate(min_year=min(year)) %>%
  filter(year==min_year) %>%
  xtabs(~year, data=.)

# Distribution of "year of POS test"
dat_prc %>%
  group_by(id) %>%
  filter(!is.na(ResultDate) & HIVResult=="P") %>%
  xtabs(~year, data=.)

# Distribution of "year of first POS test"
dat_prc %>%
  group_by(id) %>%
  filter(!is.na(ResultDate) & HIVResult=="P") %>%
  mutate(min_year=min(year)) %>%
  filter(year==min_year) %>%
  xtabs(~year, data=.)

# Out of "year at entry" years, which years had HIV tests?
dat_prc %>%
  group_by(id) %>%
  mutate(min_year=min(year)) %>%
  filter(year==min_year) %>%
  filter(!is.na(ResultDate)) %>%
  xtabs(~year, data=.)


# !!!!! New DQA: checking "initial sero" model
dat_2 <- dat %>%
  group_by(id) %>%
  mutate(min_time=min(t_end), max_time=max(t_end))
dat_3 <- dplyr::filter(dat_2, t_end==min_time)
dat_4 <- dplyr::filter(dat_3, v==1)
dat_4 %<>% dplyr::mutate(
  age_below_50 = In(age<=50),
  year_00_07 = In(t_end<=7),
  year_08_16 = In(t_end>8 & t_end<=16),
  year_17_22 =  In(t_end>16)
)

# % Positive (out of testers)
mean(dplyr::filter(dat_4, T)$u) # Overall 34%
mean(dplyr::filter(dat_4, w_1==1)$u) # Male: 24%
mean(dplyr::filter(dat_4, w_1==0)$u) # Female: 38%
length(dplyr::filter(dat_4, year_00_07==1)$u) # Year 00-07: xx%
length(dplyr::filter(dat_4, year_08_16==1)$u) # Year 08-16: xx%
length(dplyr::filter(dat_4, year_17_22==1)$u) # Year 17-22: xx%

nrow(dplyr::filter(dat_4, w_1==0 & age_below_50==1 & year_00_07==1))
nrow(dplyr::filter(dat_4, w_1==1 & age_below_50==1 & year_00_07==1))
nrow(dplyr::filter(dat_4, w_1==0 & age_below_50==1 & year_08_16==1))
nrow(dplyr::filter(dat_4, w_1==1 & age_below_50==1 & year_08_16==1))
nrow(dplyr::filter(dat_4, w_1==0 & age_below_50==1 & year_17_22==1))
nrow(dplyr::filter(dat_4, w_1==1 & age_below_50==1 & year_17_22==1))

chk(2, "DQA: END")
