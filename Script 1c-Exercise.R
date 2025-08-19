#Practical Exercise

#1 Checking Cholesterol level using the "if"

cholesterol <- 230
if(cholesterol > 240){
  print("High Cholesterol")
}

#2 Check Blood pressure status using "if" and "else"

systolic_bp <- 130

if(systolic_bp < 120){
  print("Blood Pressure is Normal")
}else{
  print("Blood Pressure is High")
}

#3 Automating Data type With for loop
Patient_info <- read.csv(file.choose())
str(Patient_info)

#Creating a copy
copy_patient <- Patient_info
View(copy_patient)

factor_cols <- copy_patient[c("gender","diagnosis", "smoker")]
str(factor_cols)

for (i in 1:ncol(factor_cols)) {
  if(is.character(factor_cols[[i]])){
    factor_cols[[i]] <- as.factor(factor_cols[[i]])
  }
}
str(factor_cols)

#converting factors to numbers codes
binary_cols <- c(factor_cols$smoker)

for (col in binary_cols) {
  if(is.factor(binary_cols[[i]])){
    binary_cols <- ifelse(binary_cols == "Yes", 1,0)
  }
}
str(binary_cols)














