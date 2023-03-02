function [equal_var] = F_test(sigma1,sigma2,n1,n2)
%Test if two measurements, given their standard devition and number of
%measurements, have equal variance. Outputs a true value if variances are
%equal.

  %F-test
    %Equal Variance?
    df_1 = n1-1;
    df_2 = n2-1;

    %X = finv(0.975,num_df,denom_df);
    p = 0.05;
    use_alpha = 1-(p/2);
    if sigma1>sigma2
        F_calc = (sigma1^2)/(sigma2^2);
        F_crit = finv(use_alpha,df_1,df_2);
    else
        F_calc = (sigma2^2)/(sigma1^2);
        F_crit = finv(use_alpha,df_2,df_1);
    end

    if F_crit<F_calc
        equal_var = false;
    else %F_crit>F_calc
        equal_var = true;
    end

end