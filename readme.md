# IHYDE

The source code of Data Driven Discovery of Cyber Physical Systems. 

## Test platform

* **Win10**, ```matlab 2017a```

## Tips
For Matlab 2018a and later version, you should set the algorithm to trust-region when you try to identify the transition logic using slr.
It can be realized by replacing the code in line 103 of function "slr_learning_l1" as:

    option = optimset('Gradobj','on','Hessian','on',...
       'MaxIter', WMaxIter, 'Display', WDisplay,'Algorithm','trust-region');

CVX toolbox (http://cvxr.com/cvx/) is need to run the code.


## Citation
Yuan, Y., Tang, X., Zhou, W. et al. Data driven discovery of cyber physical systems. Nat Commun 10, 4894 (2019).
https://www.nature.com/articles/s41467-019-12490-1
