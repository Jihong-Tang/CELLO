from sklearn.cluster import KMeans
import numpy as np

def do_clustering(df_mut_num):
    """
    """
    tmp_sum = df_mut_num.Primary + df_mut_num.Recurrent + df_mut_num.Common
    mt_cluster = np.array([df_mut_num.Primary / tmp_sum,
                        df_mut_num.Recurrent / tmp_sum, df_mut_num.Common / tmp_sum]).T
    kmeans = KMeans(n_clusters=3, max_iter=1000).fit(mt_cluster)
    df_cluster = pd.DataFrame(
        mt_cluster, columns=['Primary', 'Recurrence', 'Common'])
    clusters = np.array(kmeans.labels_)

    # set the cluster feature based on the P R C pattern
    bool_P = (clusters == clusters[df_cluster.Primary == max(df_cluster.Primary)])
    bool_R = (
        clusters == clusters[df_cluster.Recurrence == max(df_cluster.Recurrence)])
    bool_C = (clusters == clusters[df_cluster.Common == max(df_cluster.Common)])

    clusters[bool_P] = [0] * sum(bool_P)
    clusters[bool_R] = [1] * sum(bool_R)
    clusters[bool_C] = [2] * sum(bool_C)
    df_cluster = df_cluster.assign(
        Cluster=clusters).assign(id=range(len(df_cluster)))
    return(df_cluster)

