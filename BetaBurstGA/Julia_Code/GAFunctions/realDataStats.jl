# the form of this matrix is:
#    __________________________________________________________
#   | mean Duration     |   mean power     | mean number       |
#   |___________________|__________________|___________________|
#   | std Duration      | std power        | std number        |
#   |___________________|__________________|___________________|
#   | skew Duration     | skew power       | skew number       |
#   |___________________|__________________|___________________|
#   | kurtosis Duration | kurtosis power   | kurtosis number   |
#   |___________________|__________________|___________________|
#
Burst_Dur_mean = 265.61804487179484 ; Burst_Power_mean = 0.7310285446153846;
Burst_Number_mean = 1.132948441025641 ; Burst_Dur_std = 41.26250201528616;
Burst_Power_std = 0.09068969468335884; Burst_Number_std = 0.2124833051350967;
Burst_Dur_skew = -0.023487537015022993; Burst_Power_skew = 0.5331688354498735;
Burst_Number_skew = 0.7886700989825924; Burst_Dur_kurtosis = 10.08795890926448;
Burst_Power_kurtosis = -0.1976111936139313; Burst_Number_kurtosis = 2.794143075139668;
realdata_stats = [Burst_Dur_mean Burst_Power_mean Burst_Number_mean ;
                    Burst_Dur_std Burst_Power_std Burst_Number_std;
                    Burst_Dur_skew Burst_Power_skew Burst_Number_skew;
                    Burst_Dur_kurtosis Burst_Power_kurtosis Burst_Number_kurtosis;
                    ]
