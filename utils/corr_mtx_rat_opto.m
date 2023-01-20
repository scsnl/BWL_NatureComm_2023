function[sig_t_mtx_ONvsOFF, sig_t_mtx_ONvsTrans, sig_t_mtx_TransvsOFF]=corr_mtx_rat_opto(sig_t_mtx_ONvsOFF, sig_t_mtx_ONvsTrans, sig_t_mtx_TransvsOFF)

sig_t_mtx_ONvsOFF(4,[5,6])=0;
sig_t_mtx_ONvsOFF(5,[4,6])=0;
sig_t_mtx_ONvsOFF(6,[4,5,7])=0;
sig_t_mtx_ONvsOFF(7,6)=0;
sig_t_mtx_ONvsOFF(7,8)=0;
sig_t_mtx_ONvsOFF(8,7)=0;
sig_t_mtx_ONvsOFF(8,9)=0;
sig_t_mtx_ONvsOFF(9,8)=0;

sig_t_mtx_ONvsTrans(4,[5,6])=0;
sig_t_mtx_ONvsTrans(5,[4,6])=0;
sig_t_mtx_ONvsTrans(6,[4,5,7])=0;
sig_t_mtx_ONvsTrans(7,6)=0;
sig_t_mtx_ONvsTrans(7,8)=0;
sig_t_mtx_ONvsTrans(8,7)=0;
sig_t_mtx_ONvsTrans(8,9)=0;
sig_t_mtx_ONvsTrans(9,8)=0;
sig_t_mtx_ONvsTrans(5,7)=0;
sig_t_mtx_ONvsTrans(7,5)=0;

sig_t_mtx_TransvsOFF(4,[5,6])=0;
sig_t_mtx_TransvsOFF(5,[4,6])=0;
sig_t_mtx_TransvsOFF(6,[4,5,7])=0;
sig_t_mtx_TransvsOFF(7,6)=0;
sig_t_mtx_TransvsOFF(7,8)=0;
sig_t_mtx_TransvsOFF(8,7)=0;
sig_t_mtx_TransvsOFF(8,9)=0;
sig_t_mtx_TransvsOFF(9,8)=0;
sig_t_mtx_TransvsOFF(6,9)=0;
sig_t_mtx_TransvsOFF(9,6)=0;