1. Run donor/acceptor grid-search on training dataset (*Blade1* genus for *Ascomycota*)
    * This creates model for each combination of parameters `degree`, `C` and `window-size`
    
2. Validate the models on test set

3. Using the script `evaluate-models.sh` with a minor alterations, print recall&precision of each validated model to a separate files:
    * ./gridsearch-acceptor-precisions.txt
    * ./gridsearch-acceptor-recalls.txt
    * command to extract metrics from grid search using *Blade1*  data:
    
           bash evaluate_models.sh -f ../gridsearch/grid_validation_results/ascomycota/donor/ -r 0.00254 -t 0.188 -a 0.5 -g
      where parameters for precision adjustment were taken as from `./adjusting-precision-methodology.md`
      13210 true introns 3 449 328 acceptor site candidates and 5 248 934 donor site candidates makes the ratio 0.00377 and 0.00254 respectively
      0.188 is a ratio taken from TP, TN, FP and FN in result file of validation
           
4. Explore the metrics in the R script
    * The statistics and conclusions are included in the Rmd output
    * **NOTE:** take care of what metric belongs to which model (in file they are just printed on newline without structure)
    The output of the metrics should be in such format:
    
    ```    
    Metrics for model /home/anhvu/PycharmProjects/mycointrons/gridsearch/grid_validation_results/ascomycota/donor//train-splice-site-shuffeled_train_Blade1-donor-windows.csv-100-15-70-60-14066995.arien-pro.ics.muni.cz-results.txt: TPs: 11281.0, FPs: 1314.0, FNs = 1707.0
    0.895673
    Metrics for model /home/anhvu/PycharmProjects/mycointrons/gridsearch/grid_validation_results/ascomycota/donor//train-splice-site-shuffeled_train_Blade1-donor-windows.csv-100-15-70-70-14066996.arien-pro.ics.muni.cz-results.txt: TPs: 11305.0, FPs: 1278.0, FNs = 1683.0
    0.898434
    Metrics for model /home/anhvu/PycharmProjects/mycointrons/gridsearch/grid_validation_results/ascomycota/donor//train-splice-site-shuffeled_train_Blade1-donor-windows.csv-100-15-70-80-14066997.arien-pro.ics.muni.cz-results.txt: TPs: 11300.0, FPs: 1255.0, FNs = 1688.0
    0.90004
    Metrics for model /home/anhvu/PycharmProjects/mycointrons/gridsearch/grid_validation_results/ascomycota/donor//train-splice-site-shuffeled_train_Blade1-donor-windows.csv-100-20-70-60-14067004.arien-pro.ics.muni.cz-results.txt: TPs: 11261.0, FPs: 1228.0, FNs = 1727.0
    0.901673
    Metrics for model /home/anhvu/PycharmProjects/mycointrons/gridsearch/grid_validation_results/ascomycota/donor//train-splice-site-shuffeled_train_Blade1-donor-windows.csv-100-20-70-70-14067005.arien-pro.ics.muni.cz-results.txt: TPs: 11276.0, FPs: 1201.0, FN
    ```
    therefore the grid search table should have dimension `C x d x win_in`
    
    * Results can be used to determine the best parametrization for splice site models