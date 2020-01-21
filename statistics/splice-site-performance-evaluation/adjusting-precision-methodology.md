Donor/Acceptor model performance metrics cannot be used as they are - specifically precision needs to be adjusted.
This is due to the fact we trained (and tested) the models on data, that have different ratio of +/- classes, than
there is in actual data (meaning the test and train sets have more positive donor and acceptor sites, than there is
in mushrooms). To address this positive-negative examples imbalance we need to find the true ratio of positive vs. negative
AG/GT dimers in mushrooms. This file documents just that (process is demonstrated on *Ascomycota* division).

#### ACCEPTOR SITE ####
First we determine the +/- ratio in our test (or validation) set:

```shell script
>> grep -c ";+1" < ../data/ascomycota/test/shuff_aggreg_acceptor_site_test.csv 
198999 positive candidates

>> grep -c ";-1" < ../data/ascomycota/test/shuff_aggreg_acceptor_site_test.csv
300951 negative candidates
```
set the
**test-set ratio to 0.66**
---------------------------------------------------------------------------------------------------------------------
Then we estimate the true positive/negative acceptor ratio from results of 2 mushrooms:

True ratio estimation taken from *Tripe1* and *Chafi1* pipeline results (*Tripe1* for caluclation demo). First calculate how many
AG pair there are:

```shell script
>> wc -l < ../pipeline/Tripe1_results/splice-site-acceptor-dataset.csv
876865
```

where the **number of true introns on positive strand is 1563** (see `../pipeline/Tripe1_results/pipeline.output`)
Thus the ratio (in case of *Tripe1* strain) is **0.001782486** (true introns / total number of acceptor sites)

The same calculation on *Chafi1* gives ratio **0.004125975**
We take the mean of the two numbers which gives us ratio of **0.003** (cca 1:300 pos:neg donor sites).

Using the numbers we then run performance script like so:

```shell script
>> bash evaluate_models.sh -f ~/PycharmProjects/mycointrons/pipeline/bestmodels/ascomycota/classify-splice-sites-acceptor-model.log -r 0.003 -t 0.66 -a 0.5
Score 0.462433, recall 0.826678, adjusted precision 0.155342
```

#### DONOR SITE 
First we determine the +/- ratio in our test set:

```shell script
>> grep -c ";+1" < ../data/ascomycota/test/shuff_aggreg_donor_site_test.csv
211187 positive candidates

>> grep -c ";-1" < ../data/ascomycota/test/shuff_aggreg_donor_site_test.csv
300951 negative candidates
```

set the **test-set ratio to 0.7017**

---------------------------------------------------------------------------------------------------------------------
True ratio estimation taken from combined *Tripe1* and *Chafi1* pipeline results (all splice sites):
```shell script
>> wc -l < ../pipeline/Tripe1_results/splice-site-donor-dataset.csv
819784
```

where **number of true introns on positive strand is 1563** (see ../pipeline/Tripe1_results/pipeline.output)
thus the true ratio estimation (in case of *Tripe1* strain) is **0.00190668** (true introns/total number of donor sites)

The same calculation on *Chafi1* gives ratio **0.00439363**
We take the mean of the two numbers which gives us ratio of **0.00315** (cca 1:300 pos:neg donor sites).

