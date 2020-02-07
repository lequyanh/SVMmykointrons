import pandas as pd

nn_output = f'nn_acceptor_output.csv'
svm_output = f'svm_acceptor_output.csv'

nn_classification = pd.read_csv(nn_output, delimiter=';')
svm_classification = pd.read_csv(svm_output, delimiter=';')

# Adjust range so the sequences of two files are equally long
nn_classification['sequence'] = nn_classification[['sequence']].applymap(lambda seq: seq[128:270])
# nn_classification.to_csv("nn_acceptor_trimmed.csv", index=False, sep=";")

nn_classification.drop_duplicates(inplace=True)
svm_classification.drop_duplicates(inplace=True)

joined_results = pd.merge(nn_classification, svm_classification, on='sequence')
positives_df = joined_results[joined_results['label_x'] == 1]

fn_x = positives_df[positives_df['pred_x'] == -1]
fn_y = positives_df[positives_df['pred_y'] == -1]

tp_x = positives_df[positives_df['pred_x'] == 1]
tp_y = positives_df[positives_df['pred_y'] == 1]

tp_tp_intersect = positives_df[(positives_df['pred_x'] == 1) & (positives_df['pred_y'] == 1)]
print(f'True positives intersection {len(tp_tp_intersect)}. '
      f'I.e {100 * len(tp_tp_intersect) / len(tp_x)}% of model X '
      f'and {100 * len(tp_tp_intersect) / len(tp_y)}% of model Y results')

fn_fn_intersect = positives_df[(positives_df['pred_x'] == -1) & (positives_df['pred_y'] == -1)]
print(f'False negatives intersection {len(fn_fn_intersect)}. '
      f'The intersection makes {100 * len(fn_fn_intersect) / len(fn_x)}% of model X false negatives '
      f'and {100 * len(fn_fn_intersect) / len(fn_y)}% of model Y false negatives')

fn_tp_intersect = positives_df[(positives_df['pred_x'] == -1) & (positives_df['pred_y'] == 1)]
print(f'FN vs TP intersection {len(fn_tp_intersect)}. '
      f'I.e {100 * len(fn_tp_intersect) / len(fn_x)}% of Xs model FN is found by model Y ')

tp_fn_intersect = positives_df[(positives_df['pred_x'] == 1) & (positives_df['pred_y'] == -1)]
print(f'FN vs TP intersection {len(tp_fn_intersect)}. '
      f'I.e {100 * len(tp_fn_intersect) / len(fn_y)}% of Ys model FN is found by model X ')

negatives_df = joined_results[joined_results['label_x'] == -1]
fp_x = negatives_df[negatives_df['pred_x'] == 1]
fp_y = negatives_df[negatives_df['pred_y'] == 1]

fp_fp_intersect = negatives_df[(negatives_df['pred_x'] == 1) & (negatives_df['pred_y'] == 1)]
print(f'False positives intersection {len(fp_fp_intersect)}. '
      f'The intersection makes {100 * len(fp_fp_intersect) / len(fp_x)}% of model X false positives '
      f'and {100 * len(fp_fp_intersect) / len(fp_y)}% of model Y false positives')
