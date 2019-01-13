import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
from sklearn.svm import LinearSVC, SVC
from sklearn.linear_model import SGDClassifier, LogisticRegression
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV, cross_val_score, StratifiedKFold
from sklearn.feature_selection import SelectFromModel, RFE, RFECV
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
import plotnine


def grid_search_models(X,y):

    # get only exons 4-12

    X2 = X[:,3:12]
    X_train, X_test, y_train, y_test = train_test_split(X2,y,test_size=0.3)

    #SVM

    svc = SVC()
    param_grid = {'C':[0.5,1,2,3,5,6,7,8,9,10],'kernel':['rbf','linear','poly','sigmoid'],'degree':[2,3,4,5,6]}
    grid_search_svc = GridSearchCV(svc, param_grid,
                               scoring='accuracy')
    grid_search_svc.fit(X_train, y_train)

    #logistic regression

    lr = LogisticRegression()
    param_grid = {'penalty':['l1','l2'],'C':[0.5,1,2,3,4,5,8,10]}
    grid_search_lr = GridSearchCV(lr, param_grid,
                               scoring='accuracy')
    grid_search_lr.fit(X_train, y_train)

    #decision tree

    dt = DecisionTreeClassifier()
    param_grid = {'max_depth': [3, 10, 20, 30], 'max_leaf_nodes': [2, 4, 6, 8],'min_samples_leaf':[1,2,3],'min_samples_split':[2,4,6]}
    grid_search_dt = RandomizedSearchCV(dt, param_grid, cv=10,
                               scoring='accuracy')
    grid_search_dt.fit(X_train, y_train)

    # plot performances

    data = {
        'Model':['SVM']*10 + ['LogisticRegression']*10 + ['DecisionTree']*10,
        'Accuracy':list(cross_val_score(grid_search_svc.best_estimator_,X_train,y_train,cv=10)) + \
        list(cross_val_score(grid_search_lr.best_estimator_,X_train,y_train,cv=10)) + \
        list(cross_val_score(grid_search_dt.best_estimator_,X_train,y_train,cv=10))
    }
    data = pd.DataFrame(data)
    data['Model'] = pd.Categorical(data['Model'], categories=['SVM','LogisticRegression','DecisionTree'], ordered=True)

    p = plotnine.ggplot(data,plotnine.aes('Model','Accuracy')) + plotnine.geom_boxplot() + plotnine.ylim(0,1)
    p.save('./plots/tumor_genotype_prediction/accuracy-model.png')


def pca_analysis(sample_info):

    X = pd.read_csv('/Users/charlesmurphy/Dropbox/Research/0914_hui/results/RNAseq/Cufflinks/genes-symbols-fpkm.csv',index_col=0)
    common = list(set(sample_info['sample_ID'].tolist()).intersection(set(X.columns.tolist())))

    X = X[list(common)]
    sample_info = sample_info.loc[common,:]

    # do the PCA

    X = X[(X.T != 0).any()]
    X = X[X.std(axis=1)>1]
    X = np.log(X + 0.1) # log transorm is necessary to avoid heteroskedasticity
    X = scale(X,axis=1)


    pca = PCA(n_components=2)
    pca.fit(X)

    genotype = sample_info['brca1'].values

    fig, (ax1, ax2) = plt.subplots(1,2)
    fig.set_size_inches(10,5)
    ax1.scatter(
        pca.components_[0,genotype=='flox/flox'],
        pca.components_[1,genotype=='flox/flox'],
        color='red',label="Brca1 DEL")
    ax1.scatter(
        pca.components_[0,genotype=='WT'],
        pca.components_[1,genotype=='WT'],
        color='black',label="Brca1 WT")
    #ax1.legend(bbox_to_anchor=(1, 0.5))
    ax1.set_xlabel('PC1 ({0}%)'.format(round(pca.explained_variance_ratio_[0]*100)))
    ax1.set_ylabel('PC2 ({0}%)'.format(round(pca.explained_variance_ratio_[1]*100)))
    ax1.set_title('All samples')

    # remove potentially mislabeled samples and format the data

    tmp = pd.DataFrame([
        list(common),
        list(genotype),
        list(pca.components_[0]),
        list(pca.components_[1])
    ]).T

    tmp_nooutliers = tmp[((tmp[1]=='flox/flox') & (tmp[2]<=0.025)) | ((tmp[1]=='WT') & (tmp[2]>=0.025))]

    ax2.scatter(
        tmp_nooutliers.loc[tmp_nooutliers[1]=='flox/flox',2],
        tmp_nooutliers.loc[tmp_nooutliers[1]=='flox/flox',3],
        color='red',label="Brca1 DEL")
    ax2.scatter(
        tmp_nooutliers.loc[tmp_nooutliers[1]=='WT',2],
        tmp_nooutliers.loc[tmp_nooutliers[1]=='WT',3],
        color='black',label="Brca1 WT")
    ax2.legend(bbox_to_anchor=(1, 0.5))
    ax2.set_xlabel("PC1 ({0}%)".format(round(pca.explained_variance_ratio_[0]*100)))
    ax2.set_ylabel("PC2 ({0}%)".format(round(pca.explained_variance_ratio_[1]*100)))
    ax2.set_title('Potential outliers removed')
    plt.tight_layout()
    plt.savefig('./plots/tumor_genotype_prediction/pca.png')
    plt.clf()
    plt.close()

    return tmp_nooutliers

#load sample information

sample_info = pd.read_excel('/Users/charlesmurphy/Dropbox/Research/0914_hui/data/samples.xlsx','samples')
sample_info = sample_info[(sample_info['sequencing_type']=='RNAseq') & (sample_info['tissue']!='normal breast')]
sample_info.index = sample_info['sample_ID']

#PCA analysis

tmp_nooutliers = pca_analysis(sample_info)

# plot discriminatory power of different exons

exon_counts = pd.read_csv('./data/tumor_genotype_prediction/method1/brca1.exon_counts.rnaseq.csv')
common = list(set(sample_info['sample_ID'].tolist()).intersection(set(exon_counts.columns.tolist())))
exon_counts = exon_counts[list(common)]
sample_info = sample_info.loc[common,:]

data = []

for i in exon_counts.index:
    for sample in exon_counts.columns:

        tmp = sample_info.at[sample,'brca1']
        if tmp=='flox/flox':
            tmp="DEL"

        data.append([
            sample,
            tmp,
            'exon' + str(i+1),
            exon_counts.at[i,sample]
        ])
data = pd.DataFrame(data,columns=['sample','brca1','exon','normalized_read_count'])
data['exon'] = pd.Categorical(data['exon'], categories=['exon' + str(i+1) for i in range(23)], ordered=True)

p = plotnine.ggplot(data,plotnine.aes('exon','normalized_read_count',color='brca1')) + \
    plotnine.geom_boxplot() + \
    plotnine.theme(axis_text_x=plotnine.element_text(rotation=90, hjust=1))
p.save('./plots/tumor_genotype_prediction/exon_counts.png',width=15,height=5)



exon_counts = exon_counts.T
exon_counts = exon_counts.loc[exon_counts.index.isin(tmp_nooutliers[0].tolist()),:]
X = exon_counts.values
y = sample_info.loc[exon_counts.index,"brca1"]
y[y=='flox/flox']='1'
y[y=='WT']='0'
y = y.astype(int).values

# train the models

X_train, X_test, y_train, y_test = train_test_split(X,y,test_size=0.3)
#grid_search_models(X,y)

# search for optimal number of features

best_model = SVC(C=6,kernel='linear',degree=2)
rfecv = RFECV(estimator=best_model, step=1, cv=StratifiedKFold(10),
              scoring='accuracy')
rfecv.fit(X,y)

plt.figure()
plt.xlabel("Number of features selected")
plt.ylabel("Accuracy (10-fold CV)")
plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_)
plt.savefig('./plots/tumor_genotype_prediction/accuracy-num-features.png')
plt.clf()

#predict brca1 status for potential outlier samples

X = X[:,rfecv.ranking_==1]
best_model.fit(X,y)

exon_counts = pd.read_csv('./data/tumor_genotype_prediction/method1/brca1.exon_counts.rnaseq.csv')
exon_counts = exon_counts[list(common)]
exon_counts = exon_counts.T
exon_counts = exon_counts.loc[~exon_counts.index.isin(tmp_nooutliers[0].tolist()),:]
X_pred = exon_counts.values[:,rfecv.ranking_==1]

pred = best_model.predict(X_pred)
pred = ['flox/flox' if i==1 else 'WT' for i in pred]
listed = sample_info.loc[exon_counts.index.tolist(),"brca1"].tolist()
list(zip(listed,pred))
