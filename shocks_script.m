% Shocks med NQ

mednq = [med_NQ_2;
         med_NQ_3; 
         med_NQ_4; 
         med_NQ_5; 
         med_NQ_6; 
         med_NQ_7; 
         med_NQ_8; 
         med_NQ_9; 
         med_NQ_10; 
         med_NQ_11; 
         med_NQ_12; 
         med_NQ_13; 
         med_NQ_14; 
         med_NQ_15; 
         med_NQ_16; 
         med_NQ_17];

% shocks med Q
medq = [med_Q_2; 
        med_Q_3; 
        med_Q_4; 
        med_Q_5; 
        med_Q_6; 
        med_Q_7; 
        med_Q_8; 
        med_Q_9; 
        med_Q_10; 
        med_Q_11; 
        med_Q_12; 
        med_Q_13; 
        med_Q_14; 
        med_Q_15; 
        med_Q_16; 
        med_Q_17];


% shocks migration NQ
mignq = [mig_NQ_2; 
         mig_NQ_3; 
         mig_NQ_4; 
         mig_NQ_5; 
         mig_NQ_6; 
         mig_NQ_7; 
         mig_NQ_8; 
         mig_NQ_9; 
         mig_NQ_10; 
         mig_NQ_11; 
         mig_NQ_12; 
         mig_NQ_13; 
         mig_NQ_14; 
         mig_NQ_15; 
         mig_NQ_16; 
         mig_NQ_17];


% shocks migration Q
migq = [mig_Q_2; 
        mig_Q_3; 
        mig_Q_4; 
        mig_Q_5; 
        mig_Q_6; 
        mig_Q_7; 
        mig_Q_8; 
        mig_Q_9; 
        mig_Q_10; 
        mig_Q_11; 
        mig_Q_12; 
        mig_Q_13; 
        mig_Q_14; 
        mig_Q_15; 
        mig_Q_16; 
        mig_Q_17];



a_nq = [a_NQ_1;
        a_NQ_2;
        a_NQ_3;
        a_NQ_4;
        a_NQ_5;
        a_NQ_6;
        a_NQ_7;
        a_NQ_8;
        a_NQ_9];

a_q = [a_Q_1;
       a_Q_2;
       a_Q_3;
       a_Q_4;
       a_Q_5;
       a_Q_6;
       a_Q_7;
       a_Q_8;
       a_Q_9];

save('mig_q_shocks.mat', 'migq')
save('mig_nq_shocks.mat', 'mignq')
save('health_q_shocks.mat', 'medq')
save('health_nq_shocks.mat', 'mednq')


plot(1900+(1:40)*5, mednq)
plot(1900+(1:40)*5, medq)
plot(1900+(1:40)*5, mignq)
plot(1900+(1:40)*5, migq)