##############################
############ exec ############
##############################

source('xml_function.R')

#### data generate ####
## train_data
set.seed(1)
n = 100L
input_dim = 5L
output_dim = 3L
x_data = matrix(rnorm(n * input_dim), ncol = input_dim)
beta_star = matrix(rnorm(input_dim * output_dim), nrow = input_dim, ncol = output_dim)
y_prop =  exp(x_data %*% beta_star)/(1 + exp(x_data %*% beta_star))
y_data = matrix(rbinom(n * output_dim, size = 1, prob = as.numeric(y_prop)), ncol = output_dim)

## test_data
test_n = 10L
x_test_data = matrix(rnorm(test_n * input_dim), ncol = input_dim)
y_test_prop =  exp(x_test_data %*% beta_star)/(1 + exp(x_test_data %*% beta_star))
y_test_data = matrix(rbinom(test_n * output_dim, size = 1, prob = as.numeric(y_test_prop)), ncol = output_dim)


#### code execution ####
tree_list = fastXML_func(x = x_data, y = y_data, K = 1, MaxLeaf = 5L)
predict_list = fastXML_predict(tree_list = tree_list, x = x_data)
predict_test_list = fastXML_predict(tree_list = tree_list, x = x_test_data)



