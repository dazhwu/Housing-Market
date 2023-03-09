rm(list = ls(all = TRUE))

library(XML)
library(data.table)
library(readxl)
library(stringr)
library(fedmatch)
library(ivx)

library(NetworkToolbox)
library(stargazer)
#library(Hmisc)
library(ggplot2)
library(reshape2)
library(erer)
library(tseries)
library(frequencyConnectedness)
library(Spillover)
library(fastDummies)

library(HDGCvar)
library(igraph)



DATE_LEVEL=4  #quarterly
source("../../Common_Lib/Common_Functions.R")
source("FUN_multivariate.R")

check_gap<- function(dt, id, seq){
  temp=dt[,.N, by=id]
  temp2=dt[, .(num_rec=.N, max_seq=max(seq), min_seq=min(seq)), by=id]
  temp2[, gap:=ifelse(max_seq-min_seq+1==num_rec, TRUE, FALSE)]
  return (temp2)


}

# https://www.census.gov/manufacturing/m3/historical_data/index.html
format_name <- function(dt, num_char) {
  dt[, area_name := str_replace_all(area_name, " ", "")]
  dt[, area_name := str_replace_all(area_name, ",", "")]
  dt[, area_name := substring(area_name, 1, num_char)]

}

get_match <- function(dt, maxD=0.02) {
  num_char = 20
  format_name(dt, num_char)

  area = fread("./price_index/cu_area.txt")
  price_index_area = area[area_code %like% "S[0-9][0-9][A-E]", .(area_code,
                                                                 area_name)]
  format_name(price_index_area, num_char)
  
  fuzzy_result <- merge_plus(data2 = dt,
                             data1 = price_index_area,
                             by.x = "area_name",
                             by.y = "area_name", match_type = "fuzzy",
                             fuzzy_settings = build_fuzzy_settings(maxDist = maxD),
                             unique_key_2 = "id",
                             unique_key_1 = "area_code")
  return(fuzzy_result$matches)

}

download_regional_employment <- function() {
  #https://download.bls.gov/pub/time.series/sm/sm.data.1.AllData
  #Regional_employment
  download.file("https://download.bls.gov/pub/time.series/sm/sm.data.1.AllData", "./regional_employment/regional_employment.txt", mode = "wb")
  download.file("https://download.bls.gov/pub/time.series/sm/sm.area", "./regional_employment/area.txt", mode = "wb")
}
process_regional_employment <- function(){
  reg_emp_area = fread("./regional_employment/area.txt", colClasses = c("character", "character"))
  reg_emp_area = reg_emp_area[, .(area_code, area_name)]
  setnames(reg_emp_area, "area_code", "id")

  mat = get_match(reg_emp_area)[, .(id, area_code)]

  reg_emp = fread("./regional_employment/regional_employment.txt")
  reg_emp[, id := substring(series_id, 6, 10)]
  reg_emp[, data_type := substring(series_id, 19, 20)]
  reg_emp[, Month := as.integer(substring(period, 2, 3))]
  reg_emp[, industry_code := substring(series_id, 11, 18)]
  reg_emp[, season := substring(series_id, 3, 3)]
  # data type 01	All Employees, In Thousands
  #11	Average Weekly Earnings of All Employees, In Dollars

  reg_emp = reg_emp[Month <= 12]


  MSA_emp_data = reg_emp[mat, on = c("id"), nomatch = 0]
  MSA_emp_data = MSA_emp_data[, .(area_code, year, Month, data_type, value,
                                  season, industry_code)]

  num_emp = MSA_emp_data[data_type == "01" &
                           industry_code == "00000000" &
                           season == "S"]  # (data_type=="01" |
  # data_type=="11") & industry_code=="00000000"]   #01 all employees
  num_emp[, value:=as.numeric(value)]
  setnames(num_emp, "year", "Year")
  num_emp=num_emp[,.(area_code, Year, Month, value)]

  if(DATE_LEVEL==4){
    num_emp=switch2quart(num_emp)
  }

  #hourly_earning=MSA_emp_data[data_type=="03" & industry_code=="05000000"]

  defusion = MSA_emp_data[data_type %like% "2*" &
                            industry_code == "00000000" &
                            season == "S"]

  setnames(num_emp, "value", "emp")

  fwrite(num_emp, "emp.csv")


}


download_price_index <- function() {
  download.file("https://download.bls.gov/pub/time.series/cu/cu.area", "./price_index/cu_area.txt", mode = "wb")
  download.file("https://download.bls.gov/pub/time.series/cu/cu.data.0.Current", "./price_index/cu_current.txt", mode = "wb")
  #download.file("https://download.bls.gov/pub/time.series/cu/cu.data.1
  # .AllItems", "./price_index/cu_all_items.txt", mode="wb")
  price_data = fread("./price_index/cu_current.txt")
  price_data[, area_code := substring(series_id, 5, 8)]
  price_data[, Month := as.numeric(substring(period, 2, 3))]
  setnames(price_data, "year","Year")
  price_data[, p := substring(series_id, 4, 4)]  #R	Monthly     S	Semi-Annual
  fwrite(price_data, "./price_index/cu_current.txt")
}

process_price_index <- function(item_code, item_name){
  price_data = fread("./price_index/cu_current.txt")
  price_data=price_data[area_code %like% "S[0-9][0-9][A-E]"]
  price_data=price_data[series_id %like% paste0(".*", item_code, ".*")]

  price_data_monthly=price_data[p=="R" & Month <=12][, .(area_code, Year, Month, value)]
  price_data_semi_annual=price_data[p=="S" & Month <=2][, .(area_code, Year, Month, value)]


  #price_data_rent_primary = price_data_monthly[series_id %like% ".*SEHA.*", .(area_code, Year, Month, value)]

  price_data_Q=switch2quart(price_data_monthly)
  setnames(price_data_Q, "value", item_name)

  temp_1=price_data_semi_annual[, Quarter:=(Month-1)*2+1]
  temp_2=copy(price_data_semi_annual)[, Quarter:=(Month)*2]
  price_data_semi_Q = funion(temp_1, temp_2)


  #setnames(price_data_rent_semi_annual, "value", "rent")

  temp=price_data_Q[price_data_semi_Q, on=c("area_code", "Year", "Quarter")]
  new_data=temp[is.na(get(item_name))][,.(area_code, Year, Quarter, value)]
  setnames(new_data, "value",item_name)

  tbr=funion(price_data_Q, new_data)

  return (tbr)


}

# process_price_index <- function(){
#
#
#
#
#   #SEHA	Rent of primary residence
#
#   #SEHC	Owners' equivalent rent of residences	2	T	142
#   #SEHC01	Owners' equivalent rent of primary residence
#
#   #SA0L2	All items less shelter
#
#   #price_data=price_data[area_code !="0000"]
#   price_data=price_data[area_code %like% "S[0-9][0-9][A-E]"]
#
#   price_data_monthly=price_data[p=="R" & Month <=12]
#   price_data_semi_annual=price_data[p=="S" & Month <=2]
#
#   #price_data_rent_primary = price_data_monthly[series_id %like% ".*SEHA.*", .(area_code, Year, Month, value)]
#
#   price_data_rent = price_data_monthly[series_id %like% ".*SEHA.*", .(area_code, Year, Month, value)]
#
#   price_data_rent_Q=switch2quart(price_data_rent)
#   setnames(price_data_rent_Q, "value", "rent")
#
#   price_data_rent_semi= price_data_semi_annual[series_id %like% ".*SEHA.*", .(area_code, Year, Month, value)]
#   temp_1=price_data_rent_semi[, Quarter:=(Month-1)*2+1]
#   temp_2=copy(price_data_rent_semi)[, Quarter:=(Month)*2]
#   price_data_rent_semi_Q = funion(temp_1, temp_2)[,.(area_code, Year, Quarter, value)]
#
#
#   #setnames(price_data_rent_semi_annual, "value", "rent")
#
#   temp=price_data_rent_Q[price_data_rent_semi_Q, on=c("area_code", "Year", "Quarter")]
#   new_data=temp[is.na(rent)][,.(area_code, Year, Quarter, value)]
#   setnames(new_data, "value","rent")
#
#   price_data_rent_Q=funion(price_data_rent_Q, new_data)
#
#   fwrite(price_data_rent_Q, "rent.csv")
#
#   # price_data_history = price_data_rent_primary[, .(st = min(Year), end_y =max(Year)), by = .(area_code)]
#   # price_data_history = price_data_history[area, on = .(area_code)]
#   # price_data_history = price_data_history[st <= 2000 &
#   #                                           end_y == 2022 &
#   #                                           display_level == 1 &
#   #                                           area_code %like% "...[a-zA-Z]"]
#
#
#
#   price_data_CPI = price_data_monthly[series_id %like% ".*SA0L2.*", .(area_code,
#                                                               Year, Month,
#                                                               value)]
#   price_data_CPI_Q=switch2quart(price_data_CPI)
#   setnames(price_data_CPI_Q, "value", "cpi")
#
#
#   price_data_CPI_semi= price_data_semi_annual[series_id %like% ".*SEHA.*", .(area_code, Year, Month, value)]
#   temp_1=price_data_CPI_semi[, Quarter:=(Month-1)*2+1]
#   temp_2=copy(price_data_CPI_semi)[, Quarter:=(Month)*2]
#   price_data_CPI_semi_Q = funion(temp_1, temp_2)[,.(area_code, Year, Quarter, value)]
#
#
#   #setnames(price_data_rent_semi_annual, "value", "rent")
#
#   temp=price_data_CPI_Q[price_data_CPI_semi_Q, on=c("area_code", "Year", "Quarter")]
#   new_data=temp[is.na(cpi)][,.(area_code, Year, Quarter, value)]
#   setnames(new_data, "value","cpi")
#
#   price_data_CPI_Q=funion(price_data_CPI_Q, new_data)
#
#
#
#   fwrite(price_data_CPI_Q, "cpi.csv")
#
#
#
#
# }


switch2quart<-function(dt){
  dt=dt[Month<=12 & Month >=1]
  dt[, Quarter:=as.integer((Month-1)/3+1)]
  dt=dt[, .(value=mean(value)), by=.(area_code, Year, Quarter)]
  return(dt)


}

download_unemployment <- function() {

  download.file("https://download.bls.gov/pub/time.series/la/la.series", "./Unemployment/la_seris.txt", mode = "wb")
  #https://download.bls.gov/pub/time.series/la/la.area
  download.file("https://download.bls.gov/pub/time.series/la/la.area", "./Unemployment/la_area.txt", mode = "wb")
  #https://download.bls.gov/pub/time.series/la/la.area_type
  download.file("https://download.bls.gov/pub/time.series/la/la.area_type", "./Unemployment/la_area_type.txt", mode = "wb")
  download.file("https://download.bls.gov/pub/time.series/la/la.data.1.CurrentS", "./Unemployment/Season_Adj.txt", mode = "wb")
  #https://download.bls.gov/pub/time.series/la/la.data.5.RegionDivisionS

  download.file("https://download.bls.gov/pub/time.series/la/la.data.4.RegionDivisionU", "./Unemployment/RegionU.txt", mode = "wb")
  download.file("https://download.bls.gov/pub/time.series/la/la.data.5.RegionDivisionS", "./Unemployment/RegionS.txt", mode = "wb")

  data_files = c("00-04", "05-09", "10-14", "15-19", "20-24", "90-94", "95-99")

  index = 1
  for (f in data_files) {
    source = paste0("https://download.bls.gov/pub/time.series/la/la.data.0.CurrentU", f)
    dest = paste0("./Unemployment/", f, ".txt")
    download.file(source, dest, mode = "wb")
    if (index == 1) {
      d_f = fread(dest)[,.(series_id, year, period, value)]
    }else {
      temp = fread(dest)[,.(series_id, year, period, value)]
      temp[,value:=as.numeric(value)]
      d_f = funion(d_f, temp)

    }
    index=index+1

  }
  fwrite(d_f, "./Unemployment/unemp_main.csv")
}

process_unemployment<-function(measureCode, item_name){
  d_f=fread("./Unemployment/unemp_main.csv")
  series = fread("./Unemployment/la_seris.txt")
  #measure_code	measure_text
  # 03	unemployment rate
  # 04	unemployment
  # 05	employment
  # 06	labor force
  # 07	employment-population ratio
  # 08	labor force participation rate
  # 09	civilian noninstitutional population

  series = series[measure_code == measureCode]

  area = fread("./Unemployment/la_area.txt")
  area = area[, .(area_code, area_text)]
  setnames(area, c("area_code", "area_text"), c("id", "area_name"))
  mat = get_match(area)[, .(id, area_code)]


  setnames(series, "area_code", "id")
  series = series[mat, on = "id"]
  d_f1 = d_f[series, on = "series_id", nomatch = 0]
  d_f1[, Month := as.integer(substring(period, 2, 3))]


  setnames(d_f1, "year","Year")
  d_f1 = d_f1[, .(area_code, Year, Month, value)]
  if (DATE_LEVEL==4){
    d_f1=switch2quart(d_f1)
  }

  setnames(d_f1, "value", item_name)


  return (d_f1)
}

fmac_housing_index<-function(){
  #https://www.freddiemac.com/fmac-resources/research/docs/fmhpi_master_file.csv
  download.file("https://www.freddiemac.com/fmac-resources/research/docs/fmhpi_master_file.csv", "./housing_index/fmac_housing_index.csv", mode = "wb")
  fmac=fread("./housing_index/fmac_housing_index.csv")
  setnames(fmac,c("GEO_Code", "GEO_Name"), c("id", "area_name"))
  area=unique(fmac[GEO_Type=="CBSA",.(id, area_name)])
  mat=get_match(area)[,.(id, area_code)]
  fmac=fmac[mat, on="id", nomatch=0]

  fwrite(fmac[,.(area_code, Year, Month, Index_NSA, Index_SA)], "fmac_index.csv")

}

process_permits <- function() {

  permit = fread("housing_permits.txt")
  #num_char=18
  permit[, id := seq(1, nrow(permit))]
  setnames(permit, "Area_Name", "area_name")


  index = 1

  for (y in 2004:2022) {
    for (m in 1:12) {
      if ((y < 2022) | (y == 2022 & m <= 9)) {
        permit_area = permit[Year == y & Month == m, .(id, area_name)]
        mat = get_match(permit_area,0.25)[, .(id, area_code)]

        match3 = permit[mat, on = "id", nomatch = 0]
        match3 = match3[, .(area_code, Year, Month, value)]

        if (index == 1)
          combined = match3
        else
          combined = funion(combined, match3)
        index = index + 1
      }


    }
  }
  
  combined[,seq:=Month+(Year-1980)*12]
  setkeyv(combined, c("area_code","seq"))
  combined[, value:=frollmean(value, n=12), by=.(area_code)]
  if (DATE_LEVEL==4){
          combined=switch2quart(combined)
  }

  setnames(combined, "value","num_permits")
  fwrite(combined, "temp_permit.csv")

  agg = combined[, .N, by = .(area_code, Year)]
}

process_quarterly_housing_index<-function(){
  housing_index=as.data.table(read_excel("HPI_EXP_metro.xls"))
  setnames(housing_index,c("city","metro_name"), c("id","area_name"))
  housing_index[,area_name:=str_replace_all(area_name,"(MSAD)","")]
  housing_area=unique(housing_index[,.(id, area_name)])
  mat=get_match(housing_area,0.25)[,.(id, area_code)]
  housing_index=housing_index[mat, on="id", nomatch=0]
  setnames(housing_index, c("yr", "qtr"), c("Year","Quarter"))
  fwrite(housing_index[,.(area_code, Year, Quarter, index_nsa, index_sa)], "housing_index.csv")
}

process_SP500 <-function(){
  sp500=as.data.table(read_xlsx("./Macro_Fundamental/SP500.xlsx"))
  sp500[, Year:=as.integer(format(Date,'%Y'))]
  sp500[, Month:=as.integer(format(Date,'%m'))]
  sp500[, Quarter:=as.integer((Month-1)/3+1)]
  sp500=sp500[,.(Close=mean(Close)), by=.(Year, Quarter)]
  return (sp500)

}

decompose <- function(){

}

process_permits()
#process_price_index()

#process_unemployment()
process_regional_employment()
process_quarterly_housing_index()
sp500=process_SP500()


permits=fread("temp_permit.csv")
unemp=process_unemployment(4, "une")
unemp_rate=process_unemployment(3, "une_rate")


rent=process_price_index("SEHA", "rent")
cpi=process_price_index("SA0L2", "cpi")
fmac_housing=fread("fmac_index.csv")
quart_housing=fread("housing_index.csv")
emp=fread("emp.csv")

inflation_exp=fread("EXPINF1YR.csv")  #https://fred.stlouisfed.org/series/EXPINF10YR
inflation_exp[, real_date:=as.Date(DATE,format="%m/%d/%Y")]
inflation_exp[, Year:=year(real_date)+1]
inflation_exp[, Month:=month(real_date)]
inflation_exp[, Quarter:=as.integer((Month-1)/3+1)]
inflation_exp_Q=inflation_exp[, .(inflation_exp=mean(EXPINF1YR)), by=.(Year, Quarter)]

federal_fund=fread("FEDFUNDS.csv")
federal_fund[, real_date:=as.Date(DATE,format="%m/%d/%Y")]
federal_fund[, Year:=year(real_date)+1]
federal_fund[, Month:=month(real_date)]
federal_fund[, Quarter:=as.integer((Month-1)/3+1)]
federal_fund_Q=federal_fund[, .(federal_fund=mean(FEDFUNDS)), by=.(Year, Quarter)]



mortgage=as.data.table(read_excel("./Mortgage/MORTGAGE30US.xls", range="A11:B2705"))
mortgage[, observation_date:=as.character(observation_date)]
mortgage[, Year:=as.integer(substring(observation_date, 1, 4))]
mortgage[, Month:=as.integer(substring(observation_date, 6, 7))]
mortgage[, Quarter:=as.integer((Month-1)/3+1)]
if (DATE_LEVEL==4){
  mortgage=mortgage[, .(m_rate=mean(MORTGAGE30US)), by=.(Year, Quarter)]
}else{
  mortagage=mortgage[, .(m_rate=mean(MORTGAGE30US)), by=.(Year, Month)]
}

dts=list(unemp, unemp_rate, rent, cpi, quart_housing, emp, permits)

# for(dt in dts){
#   dt[, seq:=Quarter+(Year-1980)*4]
#   temp=check_gap(dt, "area_code","seq")
#   print(temp)
#   temp
# }
index=1

for (dt in dts){
  if (index==1){
    combined=dt
  }else{
    if(DATE_LEVEL==4){
      combined=combined[dt, on=.(area_code, Year, Quarter), nomatch=0]
    }else{
      combined=combined[dt, on=.(area_code, Year, Month), nomatch=0]
    }

  }

  index=index+1

}

combined=combined[mortgage, on=.(Year, Quarter), nomatch=0]
combined=combined[inflation_exp_Q, on=.(Year, Quarter), nomatch=0]
combined=combined[sp500, on=.(Year, Quarter), nomatch=0]
combined=combined[federal_fund_Q, on=.(Year, Quarter), nomatch=0]
combined[,population:=une/une_rate]
incols=colnames(combined)
incols=incols[incols!="Year"]
incols=incols[incols!="Quarter"]
incols=incols[incols!="area_code"]
outcols=paste0("log_", incols)
combined=combined[, c(outcols):=lapply(.SD, function(x){log(x)}), .SDcols=incols]
combined[, seq:=Quarter+(Year-1980)*4]
temp=check_gap(combined, "area_code","seq")
temp2=temp[min_seq<=97,.(area_code)]

combined=combined[temp2, on="area_code", nomatch=0]
combined[,.(mi=min(seq), ma=max(seq)), by=.(area_code)]
combined=combined[seq>=97]
data_table.lag(combined, c("log_index_nsa", "log_index_sa","log_m_rate","log_rent", "une_rate","federal_fund","log_cpi", "log_population", "log_num_permits", "log_emp"), "lag_",1,"area_code","seq")
combined[, y:=log_index_nsa-log_rent-(lag_log_index_nsa-lag_log_rent)]
#combined[, y:=log_index_nsa-(lag_log_index_nsa)-log_rent + lag_log_rent]
combined[, y2:=log_index_nsa-log_rent]
#combined[, y:=log_index_sa-(lag_log_index_sa)]
#combined[, y2:=log_index_sa]

combined[, d_rent:=log_rent-lag_log_rent]
combined[, int_rate:=log_m_rate-lag_log_m_rate]
combined[, d_cpi:=log_cpi - lag_log_cpi]
combined[, d_pop:=log_population - lag_log_population]
combined[, d_num_permits:=log_num_permits-lag_log_num_permits]
combined[, d_emp:=log_emp-lag_log_emp]
 combined =dummy_cols(combined, select_columns = 'Quarter')
index=1

msa_vector=as.vector(temp2$area_code)
num_msa=length(msa_vector)
integ_order=numeric(num_msa)

for (msa in msa_vector){
  print(msa)

  print(index)
  #关键在时间区间
  temp_data=combined[area_code==msa, .(y, y2, log_cpi, log_emp, log_une,  seq, log_population, log_num_permits, log_inflation_exp,federal_fund,log_m_rate, log_Close, Quarter_1, Quarter_2, Quarter_3)]
  #temp_data=combined[area_code==msa, .(y, y2, cpi, emp, une,  seq, log_population, log_num_permits, inflation_exp,federal_fund,m_rate, Quarter_1, Quarter_2, Quarter_3)]
  #temp_data=combined[area_code==msa, .(y, y2, d_rent, d_cpi, d_emp, d_pop, d_num_permits, log_m_rate)]
  temp_data=as.matrix(na.omit(temp_data))

  if (index==1){
    #
    num_ys=nrow(temp_data)
    fundamental=data.table(seq=1:(num_ys))
    bubble=data.table(seq=1:(num_ys))
  }
  if (nrow(temp_data)>=num_ys){

    x=temp_data[,3:ncol(temp_data)]
    temp_formula="y~ log_cpi + log_emp+log_population+ log_num_permits+log_m_rate"   #""
    #print(ivx_ar(as.formula(temp_formula), data=temp_data) %>% summary())
    temp_result=ivx_ar(as.formula(temp_formula), data=as.data.frame(temp_data), horizon=1, x=TRUE, y=TRUE)
    
    temp_result2=ivx_ar_fit(temp_data[,1], x, ar="forecast")
     s=temp_result$intercept  #+ temp_result$x * temp_result$coefficients
     for (i in 1:ncol(temp_result$x)){
        s=s+ temp_result$x[,i]*temp_result$coefficients[i]
      }
    # s=s-2.460754920585877
    coef=as.matrix(temp_result$coefficients)
    intercept=mean(temp_data[2:nrow(temp_data),1])-mean(as.matrix(temp_result$x) %*% coef)
    fitted=as.matrix(temp_result$x) %*% coef + intercept
    
    
    #fitted=as.matrix(temp_data$x)
    f=numeric(num_ys)
    e=numeric(num_ys)
    sum_y=0.0
    print(sum_y)
    for (i in 1:(num_ys)){
      if (i==1){
        f[i]=temp_data[1,2]
        sum_y=f[i]
      }
      else{
        #sum_y= sum_y+ temp_result$fitted[i-1]+mean(temp_data[,1])
        f[i]=f[i-1]+fitted[i]
      }
      e[i]=temp_data[i,2]-f[i]

    }

    fundamental[, eval(msa):=f]
    bubble[, eval(msa):=e]
    #print(ivx(as.formula(temp_formula), data=temp_data, horizon=1) %>% summary())
  }

  index=index+1
}

bubble=na.omit(bubble)
fundamental=na.omit(fundamental)
fwrite(fundamental, "fundamental.csv")
fwrite(bubble, "bubble.csv")

bubble=as.matrix(bubble)[2:nrow(bubble),2:ncol(bubble)]

for (i in 1:ncol(bubble)){
  print(as.vector(temp2$area_code)[i])
  vec_i=bubble[,i]
  df2=(ur.df2(vec_i, selectlags="BIC"))
  print(df2$teststat)
  ord=0
  integ_order[i]=ord
  while (df2$teststat>df2$cval[2]){
    ord=ord+1
    integ_order[i]=ord
    df2=(ur.df2(diff(vec_i,ord), selectlags="BIC"))
    print(paste0(ord, " -->",df2$teststat))
    print(df2$cval[2])

    print("----------------")
  }

  #za_test=(ur.za(bubble[,i], model="both",lag=0))
  #print(za_test@teststat)
}

  #print(ur.df2(bubble[,i], selectlags="BIC"))


VARselect(bubble[,16:17],lag.max=8,type="both")

est=VAR(bubble, p=2, type = "both")


a=spilloverDY12(est, n.ahead = 2, no.corr = F)

print(frequencyConnectedness::net(a))

b=G.spillover(est, n.ahead = 2, standardized = FALSE)
Spillover::net(b)
ts.plot(bubble)
selected_lag<-lags_upbound_BIC(bubble,p_max=10)
print(selected_lag)
network<-HDGC_VAR_all(bubble, p = selected_lag, d = 3, bound = 0.5 * nrow(bubble), 
                      parallel = TRUE, n_cores = 4)
Plot_GC_all(network, Stat_type="FS_cor",alpha=0.01, multip_corr=list(F),directed=T, layout=layout.circle, main="Network",edge.arrow.size=.2,vertex.size=5, vertex.color=c("lightblue"), vertex.frame.color="blue",vertex.label.size=2,vertex.label.color="black",vertex.label.cex=0.6, vertex.label.dist=1, edge.curved=0,cluster=list(T,5,"black",0.8,1,0)) 