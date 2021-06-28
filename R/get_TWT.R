# Get the TwoWayTalk Ls_Rr-Lr_Rs pairs
# Input:
# lr_expr_receiver_to_sender: the ligand-receptor network of interest, signaling from the receiver to the sender
# lr_expr_sender_to_receiver: the ligand-receptor network of interest, signaling from the sender to the receiver
# ligand_target_matrix_binary_receiver_to_sender: ligand_target_matrix_binary, signaling to/in sender
# ligand_target_matrix_binary_sender_to_receiver: ligand_target_matrix_binary, signaling to/in receiver
# receptor_target_matrix_binary_receiver_to_sender: receptor_target_matrix_binary, signaling to/in sender
# receptor_target_matrix_binary_sender_to_receiver: receptor_target_matrix_binary, signaling to/in receiver
# Output:
# myTWT: a list:
# myTWT$'Ls->Lr_Lr->Ls': TwoWayTalk network with Ls->Lr and Lr->Ls
# myTWT$'Rr->Lr_Rs->Ls': TwoWayTalk network with Rr->Lr and Rs->Ls

get_TWT <- function(lr_expr_receiver_to_sender, lr_expr_sender_to_receiver, 
                    ligand_target_matrix_binary_receiver_to_sender, ligand_target_matrix_binary_sender_to_receiver,
                    receptor_target_matrix_binary_receiver_to_sender, receptor_target_matrix_binary_sender_to_receiver) {
  
  TWT_p_list = list()
  for (i in 1:length(lr_expr_sender_to_receiver$eachcondition)) {
    TWT_p_list[[i]] = make_TWT_p_list(lr_expr_receiver_to_sender$eachcondition[[i]], lr_expr_sender_to_receiver$eachcondition[[i]],
                                      ligand_target_matrix_binary_receiver_to_sender, ligand_target_matrix_binary_sender_to_receiver,
                                      receptor_target_matrix_binary_receiver_to_sender, receptor_target_matrix_binary_sender_to_receiver)
  }
  names(TWT_p_list) = names(lr_expr_sender_to_receiver$eachcondition)
  
  TWT_p_list_bind_LsLrLrLs = vector()
  TWT_p_list_bind_RrLrRsLs = vector()
  for (i in 1:length(TWT_p_list)) {
    TWT_p_list_bind_LsLrLrLs = rbind(TWT_p_list_bind_LsLrLrLs, TWT_p_list[[i]]$`Ls->Lr & Lr->Ls`)
    TWT_p_list_bind_RrLrRsLs = rbind(TWT_p_list_bind_RrLrRsLs, TWT_p_list[[i]]$`Rr->Lr & Rs->Ls`)
  }
  TWT_p_list_bind_LsLrLrLs = unique(TWT_p_list_bind_LsLrLrLs)
  TWT_p_list_bind_RrLrRsLs = unique(TWT_p_list_bind_RrLrRsLs)
  colnames(TWT_p_list_bind_LsLrLrLs) = c('Lr','Rs', 'Ls', 'Rr')
  colnames(TWT_p_list_bind_RrLrRsLs) = c('Lr','Rs', 'Ls', 'Rr')
  
  myTWT = list(TWT_p_list_bind_LsLrLrLs, TWT_p_list_bind_RrLrRsLs)
  names(myTWT) = c('Ls->Lr_Lr->Ls', 'Rr->Lr_Rs->Ls')
  
  return(myTWT)
}









