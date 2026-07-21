
data(continuous_example)
data(Glass)
data(Sonar)
data(Soybean, package = "mlbench")
data(BreastCancer)
data(iris)
data(Ionosphere)
data(abalone)
data(FuelEconomy)

airquality_subset = airquality      %>% na.omit %>% slice_sample(n = max_n_dataset, replace = TRUE)
diamonds_subset = ggplot2::diamonds %>% na.omit %>% slice_sample(n = max_n_dataset, replace = TRUE) %>% mutate_if(where(is.factor), as.character)
boston_subset = MASS::Boston %>%        na.omit %>% slice_sample(n = max_n_dataset, replace = TRUE)
cars_subset = MASS::Cars93 %>%          na.omit %>% slice_sample(n = max_n_dataset, replace = TRUE) %>% select(-Make, -Model)
glass_subset = Glass %>%                na.omit %>% slice_sample(n = max_n_dataset, replace = TRUE) %>% mutate(Type = as.numeric(Type) - 1)
pima_subset = MASS::Pima.tr2 %>%        na.omit %>% slice_sample(n = max_n_dataset, replace = TRUE) %>% mutate(type = ifelse(type == "Yes", 1, 0))
sonar_subset = Sonar %>%                na.omit %>% slice_sample(n = max_n_dataset, replace = TRUE) %>% mutate(Class = ifelse(Class == "M", 1, 0))
soybean_subset = Soybean %>%            na.omit %>% slice_sample(n = max_n_dataset, replace = TRUE) %>% mutate(Class = ifelse(Class == "brown-spot", 1, 0)) %>%
	transform(plant.stand = as.numeric(plant.stand), precip = as.numeric(precip), temp = as.numeric(temp), germ = as.numeric(germ), leaf.size = as.numeric(leaf.size))
cancer_subset = BreastCancer %>%        na.omit %>% slice_sample(n = max_n_dataset, replace = TRUE) %>% select(-Id) %>% mutate(Class = ifelse(Class == "malignant", 1, 0)) %>%
	transform(Cl.thickness = as.numeric(Cl.thickness), Cell.size = factor(Cell.size, ordered = FALSE), Cell.shape = factor(Cell.shape, ordered = FALSE), Marg.adhesion = factor(Marg.adhesion, ordered = FALSE), Epith.c.size = as.numeric(Epith.c.size))
ionosphere_subset = Ionosphere %>%      na.omit %>% slice_sample(n = max_n_dataset, replace = TRUE) %>% select(-V1, -V2)
abalone_subset = abalone %>%            na.omit %>% slice_sample(n = max_n_dataset, replace = TRUE)
fuel_subset = cars2012 %>%              na.omit %>% slice_sample(n = max_n_dataset, replace = TRUE)
rm(max_n_dataset)


clamp_proportion_response = function(y){
	y = as.numeric(y)
	pmin(1, pmax(0, y))
}

finagle_different_responses_from_continuous = function(y_cont){
	y_scaled = scale(y_cont)
	list(
		continuous = y_scaled,
		incidence =  stats::plogis(as.numeric(y_scaled)),
		proportion = clamp_proportion_response((y_cont - min(y_cont) + 1e-6) / max(y_cont - min(y_cont) + 2e-6)),
		count =      round(y_cont - min(y_cont)),
		survival =   y_scaled - min(y_scaled) + 0.1,
		ordinal =    as.integer(cut(y_cont, breaks = unique(quantile(y_cont, probs = seq(0, 1, length.out = 5))), include.lowest = TRUE))
	)
}

#normalize the data
#fit a *.* model for all types
#draw y's differently in each simulation
datasets_and_response_models = list(
	pima = list(
		X = pima_subset %>% model.matrix(type ~ 0 + ., .) %>% apply(2, scale) %>% `/`(ncol(.)) %>% data.table %>% select(where(~ !any(is.na(.)))),
		y_original = list(incidence = pima_subset$type),
		beta_T = list(
			incidence =  2
		)
	),
	cancer = list(
		X = cancer_subset %>% model.matrix(Class ~ 0 + ., .) %>% apply(2, scale) %>% `/`(ncol(.)) %>% data.table %>% select(where(~ !any(is.na(.)))),
		y_original = list(incidence = cancer_subset$Class),
		beta_T = list(
			incidence =  2
		)
	),
	sonar = list(
		X = sonar_subset %>% model.matrix(Class ~ 0 + ., .) %>% apply(2, scale) %>% `/`(ncol(.)) %>% data.table %>% select(where(~ !any(is.na(.)))),
		y_original = list(incidence = sonar_subset$Class),
		beta_T = list(
			incidence =  2
		)
	),
	soybean = list(
		X = soybean_subset %>% model.matrix(Class ~ 0 + ., .) %>% apply(2, scale) %>% `/`(ncol(.)) %>% data.table %>% select(where(~ !any(is.na(.)))),
		y_original = list(incidence = soybean_subset$Class),
		beta_T = list(
			incidence =  2
		)
	),
	diamonds = list(
		X = diamonds_subset %>% model.matrix(price ~ 0 + ., .) %>% apply(2, scale) %>% `/`(ncol(.)) %>% data.table %>% select(where(~ !any(is.na(.)))),
		y_original = finagle_different_responses_from_continuous(log(diamonds_subset$price)),
		beta_T = list(
			continuous = 0.2,
			incidence =  2,
			count =      0.3
		)
	),
	pte_example = list(
		X = continuous_example$X %>% select(-treatment) %>% model.matrix(~ 0 + ., .) %>% apply(2, scale) %>% `/`(ncol(.)) %>% data.table %>% select(where(~ !any(is.na(.)))),
		y_original = finagle_different_responses_from_continuous(continuous_example$y),
		beta_T = list(
			continuous = 0.2,
			incidence =  2,
			count =      0.3,
			survival =   1
		)
	),
	ionosphere = list(
		X = ionosphere_subset %>% model.matrix(V3 ~ 0 + ., .) %>% apply(2, scale) %>% `/`(ncol(.)) %>% data.table %>% select(where(~ !any(is.na(.)))),
		y_original = finagle_different_responses_from_continuous(ionosphere_subset$V3),
		beta_T = list(
			continuous = 0.2,
			incidence =  2,
			count =      0.3,
			survival =   1
		)
	),
	abalone = list(
		X = abalone_subset %>% model.matrix(LongestShell ~ 0 + ., .) %>% apply(2, scale) %>% `/`(ncol(.)) %>% data.table %>% select(where(~ !any(is.na(.)))),
		y_original = finagle_different_responses_from_continuous(abalone_subset$LongestShell),
		beta_T = list(
			continuous = 0.2,
			incidence =  2,
			count =      0.3,
			survival =   1
		)
	),
	iris = list(
		X = iris %>% model.matrix(Sepal.Width ~ 0 + ., .) %>% apply(2, scale) %>% `/`(ncol(.)) %>% data.table %>% select(where(~ !any(is.na(.)))),
		y_original = finagle_different_responses_from_continuous(iris$Sepal.Width),
		beta_T = list(
			continuous = 0.2,
			incidence =  2,
			count =      0.3,
			survival =   1
		)
	),
	boston = list(
		X = boston_subset %>% model.matrix(medv ~ 0 + ., .) %>% apply(2, scale) %>% `/`(ncol(.)) %>% data.table %>% select(where(~ !any(is.na(.)))),
		y_original = finagle_different_responses_from_continuous(boston_subset$medv),
		beta_T = list(
			continuous = 0.2,
			incidence =  2,
			proportion = 0.3,
			count =      0.3,
			survival =   1
		)
	),
	fuel = list(
		X = fuel_subset %>% model.matrix(EngDispl ~ 0 + ., .) %>% apply(2, scale) %>% `/`(ncol(.)) %>% data.table %>% select(where(~ !any(is.na(.)))),
		y_original = finagle_different_responses_from_continuous(fuel_subset$EngDispl),
		beta_T = list(
			continuous = 0.2,
			incidence =  2,
			count =      0.3,
			survival =   1
		)
	),
	cars = list(
		X = cars_subset %>% model.matrix(Price ~ 0 + ., .) %>% apply(2, scale) %>% `/`(ncol(.)) %>% data.table %>% select(where(~ !any(is.na(.)))),
		y_original = finagle_different_responses_from_continuous(cars_subset$Price),
		beta_T = list(
			continuous = 0.2,
			incidence =  2,
			count =      0.3,
			survival =   1
		)
	),
	airquality = list(
		X = airquality_subset %>% model.matrix(Wind ~ 0 + ., .) %>% apply(2, scale) %>% `/`(ncol(.)) %>% data.table %>% select(where(~ !any(is.na(.)))),
		y_original = finagle_different_responses_from_continuous(airquality_subset$Wind),
		beta_T = list(
			continuous = 0.2,
			incidence =  2,
			proportion = 0.3,
			count =      0.3,
			survival =   1
		)
	),
	glass = list(
		X = glass_subset %>% model.matrix(Type ~ 0 + ., .) %>% apply(2, scale) %>% `/`(ncol(.)) %>% data.table %>% select(where(~ !any(is.na(.)))),
		y_original = list(count = glass_subset$Type),
		beta_T = list(
			count =      1
		)
	)
)
rm(continuous_example,Glass,Sonar,Soybean,BreastCancer,iris,Ionosphere,abalone,cars2010,cars2011,cars2012)
rm(airquality_subset,diamonds_subset,boston_subset,cars_subset,glass_subset,pima_subset,sonar_subset,soybean_subset,cancer_subset,ionosphere_subset,abalone_subset,fuel_subset)
for (dataset_name in names(datasets_and_response_models)){
	if ("proportion" %in% names(datasets_and_response_models[[dataset_name]]$y_original)) {
		datasets_and_response_models[[dataset_name]]$y_original$proportion = clamp_proportion_response(datasets_and_response_models[[dataset_name]]$y_original$proportion)
	}
}

rm(clamp_proportion_response)
rm(finagle_different_responses_from_continuous)
