#importing the patient.info csv file
data <- read.csv(file.choose())

#inspecting the data and structure of the patient.info csv file
View(data)
str(data)

#convert gender into factor
data$gender_fact <- as.factor(data$gender)
str(data)

#convert factor into numeric factor
data$gender_num <- ifelse(data$gender_fact == "Female", 1,0)
class(data$gender_num)

data$gender_num <- as.factor(data$gender_num)
class(data$gender_num)

#converting diagnosis into factor
data$diagnosis_fact <- as.factor(data$diagnosis)
str(data)


#converting smoker into a binary factor
data$smoker_fact <- as.factor(data$smoker)
data$smoker_num <- ifelse(data$smoker_fact == "Yes", 1,0)
class(data$smoker_num)

data$smoker_num <- as.factor(data$smoker_num)
class(data$smoker_num)

# save file in csv format
write.csv(data, file= "data/patient_info_clean.csv")

# Save the my workspace
save.image(file = "LanakaGaston_class_Ib_Assignment.RData")














