####################################
############# fastXML ##############
####################################

rm(list = ls())
gc(reset = T)

if(!require(dplyr)) install.packages("dplyr")
require(dplyr)

#### execution code
fastXML_func = function(x, y, K, MaxLeaf = 10L)
{
  tree_list = vector('list', K)
  predict_list = vector('list', K)
  n_id = 1:nrow(x_data)
  for(k in 1:K)
  {
    tree_list[[k]] = grow_node_recursive_func(x = x, y = y, data_inx = n_id, MaxLeaf = MaxLeaf)
  }
  return(tree_list)
}

fastXML_predict = function(tree_list, x)
{
  return(lapply(tree_list, fastXML_predict_list, x = x))
}

#### buile tree ftns
grow_node_recursive_func = function(x, y, data_inx, MaxLeaf = 10L) ## Return tree
{
  if(length(data_inx) <= MaxLeaf)
  {
    P = process_leaf_func(y, data_inx)
    w = NULL
  }else
  {
    split_node_result = split_node_func(x, y, data_inx)
    if(split_node_result$node_bug)
    {
      P = process_leaf_func(y, data_inx = split_node_result[[names(which(lapply(split_node_result[2:3], length) != 0))]])
      w = NULL
    }else{
      w = split_node_result$w
      n_left_tree = grow_node_recursive_func(x, y, split_node_result$n_left_child, MaxLeaf = MaxLeaf)
      n_right_tree = grow_node_recursive_func(x, y, split_node_result$n_right_child, MaxLeaf = MaxLeaf)
    }
  }
  if(is.null(w))
  {
    return(list(P = P))
  }else{
    return(list(w = w,
                l_tree = n_left_tree, 
                r_tree = n_right_tree))
  }
}

IL_y_func = function(s)
{
  return(1/sum(1/log(1 + 1:s)))
}

split_node_func = function(x, y, data_inx, max_iter = 1000L)
{
  set.seed(1)
  delta_next_inx = delta_inx = sample(c(-1, 1), size = length(data_inx), replace = T)
  y_id = y[data_inx, , drop = F]
  x_id = x[data_inx, , drop = F]
  d = ncol(x)
  L = ncol(y)
  
  IL_y = lapply(rowSums(y_id), IL_y_func) %>% unlist()
  IL_y_y = IL_y * y_id
  
  w_init = rep(0, d)
  # w_init = rnorm(d)
  break_stack = 0
  for(i in 1:max_iter)
  {
    r_plus_order = order(colSums(0.5 * (1 + delta_inx) * IL_y_y), decreasing = T)
    r_minus_order = order(colSums(0.5 * (1 - delta_inx) * IL_y_y), decreasing = T)
    
    v_plus = log(1 + exp(-x_id %*% w_init)) - IL_y * rowSums(y_id[, r_plus_order] / log(1 + 1:L))
    v_minus = log(1 + exp(x_id %*% w_init)) - IL_y * rowSums(y_id[, r_minus_order] / log(1 + 1:L))
    v_sign = sign(v_minus - v_plus)
    delta_next_inx[v_sign != 0] = v_sign[v_sign != 0]
    if(all(delta_next_inx == delta_inx))
    {
      w_next = optim(w_init, objective_w_func, x = x_id, delta_inx = delta_inx)$par # w update
      break_stack = break_stack + 1
      if(break_stack >= 3) break
    }else{
      w_next = w_init
    }
    delta_inx = delta_next_inx
  }
  decision_boundary = x_id %*% w_next
  n_right_child = data_inx[which(decision_boundary > 0)]
  n_left_child = data_inx[which(decision_boundary <= 0)]
  node_bug = ifelse(any(length(n_left_child) == 0, length(n_right_child) == 0), T, F)
  return(list(w = w_next, n_left_child = n_left_child, n_right_child = n_right_child, node_bug = node_bug))
}

## new idea (K-means split node)
# split_node_func = function(x, y, data_inx, max_iter = 1000L)
# {
#   set.seed(1)
#   d = ncol(x)
#   L = ncol(y)
#   
#   x_id = x_data[data_inx, ]
#   y_id = y_data[data_inx, ]
#   
#   w_init = rep(0, d)
#   break_stack = 0
#   delta_inx = F
#   for(i in 1:max_iter)
#   {
#     kmeans_clusters = kmeans(y[data_inx, ], centers = 2)
#     delta_next_inx = kmeans_clusters$cluster
#     w_init = optim(w_init, objective_w_func, x = x_id, delta_inx = delta_inx)$par # w update
#     if(all(delta_inx == delta_next_inx))
#     {
#       break_stack = break_stack + 1
#       if(break_stack >= 3) break
#     }
#     delta_inx = delta_next_inx
#   }
#   decision_boundary = x_id %*% w_init
#   n_right_child = data_inx[which(decision_boundary > 0)]
#   n_left_child = data_inx[which(decision_boundary <= 0)]
#   node_bug = ifelse(any(length(n_left_child) == 0, length(n_right_child) == 0), T, F)
#   return(list(w = w_init, n_left_child = n_left_child, n_right_child = n_right_child, node_bug = node_bug))
# }

objective_w_func = function(w, x, delta_inx)
{
  return(sum(abs(w))  + colSums(log(1 + exp(-delta_inx * x %*% w))))
}

process_leaf_func = function(y, data_inx)
{
  return(colMeans(y[data_inx, , drop = F]))
}

#### prediction ftns
fastXML_predict_list = function(tree, x)
{
  pp = apply(x, 1, fastXML_predict_each, tree = tree) %>% t()
  return(pp)
}

fastXML_predict_each = function(tree, x_vec)
{
  while(length(tree) != 1)
  {
    tree_inx = ifelse(x_vec %*% tree$w >= 0, 2, 3) %>% as.numeric()
    tree = tree[[tree_inx]]
  }
  return(tree$P)
}

