import pandas as pd
import numpy as np

def cal_HMDection(df_savi, cutoff_mut_freq, cutoff_mut_num, cutoff_HMscore):
    """
    """
    # initilize the parameters
    cutoff_varPresent = cutoff_mut_freq
    cutoff_mutLoad = cutoff_mut_num
    cutoff_HMscore = cutoff_HMscore

    case = df_savi.CaseID.drop_duplicates().tolist()

    # initilize all the useful list objects
    mutLoad_P = [0] * len(case)
    mutLoad_R = [0] * len(case)
    HMscore_P = [0] * len(case)
    HMscore_R = [0] * len(case)
    # use numpy arry for futher multiple change
    savi_P = np.array([0] * len(df_savi))
    # use numpy arry for futher multiple change
    savi_R = np.array([0] * len(df_savi))
    savi_isHM = np.array([0] * len(df_savi))

    # assign three features to gain the new HM df for further calculation
    df_savi_HM = df_savi.assign(isC2T=lambda df:
                                ((df.ref == 'C') & (df.alt == 'T')) |
                                ((df.ref == 'G') & (df.alt == 'A')))
    df_savi_HM = df_savi_HM.assign(isCC2TC=lambda df:
                                ((df.ref == 'C') & (df.alt == 'T') & (df.varSurffix == 'C')) |
                                ((df.ref == 'G') & (df.alt == 'A') & (df.varSurffix == 'G')))
    df_savi_HM = df_savi_HM.assign(isCT2TT=lambda df:
                                ((df.ref == 'C') & (df.alt == 'T') & (df.varSurffix == 'T')) |
                                ((df.ref == 'G') & (df.alt == 'A') & (df.varSurffix == 'A')))
    # change the index for further filtering work
    df_savi_HM.index = range(len(df_savi_HM))
    psedcount = 1  # parameter for HM score calculation

    for i in range(len(case)):  # for each patient
        df_patient = df_savi_HM[df_savi_HM.CaseID == case[i]]

        # primary
        df_varP = df_patient[(df_patient.Blood_freq == 0) &
                            (df_patient.Primary_freq > cutoff_varPresent)]
        num_varP = len(df_varP)
        savi_P[df_varP.index.values] = [1] * \
            len(df_varP)  # numpy array benefits here
        mutLoad_P[i] = num_varP

        # HM score calculation part IMPORTANT
        pC2T = sum(df_varP.isC2T) / (mutLoad_P[i] + psedcount)
        pCC2TC = sum(df_varP.isCC2TC) / (sum(df_varP.isC2T) + psedcount)
        pCT2TT = sum(df_varP.isCT2TT) / (sum(df_varP.isC2T) + psedcount)
        # IMPORTANT
        HMscore_P[i] = pC2T + pCC2TC + np.sign(pCC2TC - pCT2TT) * pCT2TT

        # recurrent
        df_varR = df_patient[(df_patient.Blood_freq == 0) &
                            (df_patient.Recurrent_freq > cutoff_varPresent)]
        num_varR = len(df_varR)
        savi_R[df_varR.index.values] = [1] * \
            len(df_varR)  # numpy array benefits here
        mutLoad_R[i] = num_varR

        # HM score calculation part IMPORTANT
        rC2T = sum(df_varR.isC2T) / (mutLoad_R[i] + psedcount)
        rCC2TC = sum(df_varR.isCC2TC) / (sum(df_varR.isC2T) + psedcount)
        rCT2TT = sum(df_varR.isCT2TT) / (sum(df_varR.isC2T) + psedcount)
        # IMPORTANT
        HMscore_R[i] = rC2T + rCC2TC + np.sign(rCC2TC - rCT2TT) * rCT2TT

    HMmark_P = [0] * len(case)
    HMmark_R = [0] * len(case)


    for i in range(len(case)):
        if (mutLoad_P[i] > cutoff_mutLoad) & (HMscore_P[i] > cutoff_HMscore):
            savi_isHM[df_savi_HM.CaseID == case[i]] = 1
            HMmark_P[i] = 1

        if (mutLoad_R[i] > cutoff_mutLoad) & (HMscore_R[i] > cutoff_HMscore):
            savi_isHM[df_savi_HM.CaseID == case[i]] = 1
            HMmark_R[i] = 1

    # temporary list to gain the result DataFrame
    tmp_list = []
    tmp_list.append(case * 2)  # CaseID
    tmp_arr1 = np.array(mutLoad_P) + 1
    tmp_arr2 = np.array(mutLoad_R) + 1
    tmp_list.append(np.concatenate((tmp_arr1, tmp_arr2)))  # mutation number
    tmp_list.append(HMscore_P + HMscore_R)  # HM score
    tmp_list.append(HMmark_P + HMmark_R)  # HM mark
    tmp_list.append(['P'] * len(case) + ['R'] * len(case))  # PR mark

    df_HM_Detection = pd.DataFrame(np.array(tmp_list).T, columns=[
                                'caseID', 'mutNumber', 'HMscore', 'HMmark', 'PRmark'])
    return([df_HM_Detection, df_savi_HM, savi_P, savi_R, savi_isHM])

def cal_HMFrac(df_savi_HM, savi_P, savi_R, savi_isHM):
    """
    """
    subTypes_P = df_savi_HM[savi_P == 1]
    subTypes_R = df_savi_HM[(savi_isHM == 0) & (savi_R == 1)]
    subTypes_HR = df_savi_HM[(savi_isHM == 1) & (savi_R == 1)]

    # for primary
    mut_CTGA_P = sum(((subTypes_P.ref == 'C') & (subTypes_P.alt == 'T')) |
                    ((subTypes_P.ref == 'G') & (subTypes_P.alt == 'A')))
    mut_CGGC_P = sum(((subTypes_P.ref == 'C') & (subTypes_P.alt == 'G')) |
                    ((subTypes_P.ref == 'G') & (subTypes_P.alt == 'C')))
    mut_CAGT_P = sum(((subTypes_P.ref == 'C') & (subTypes_P.alt == 'A')) |
                    ((subTypes_P.ref == 'G') & (subTypes_P.alt == 'T')))
    mut_ATTA_P = sum(((subTypes_P.ref == 'A') & (subTypes_P.alt == 'T')) |
                    ((subTypes_P.ref == 'T') & (subTypes_P.alt == 'A')))
    mut_AGTC_P = sum(((subTypes_P.ref == 'A') & (subTypes_P.alt == 'G')) |
                    ((subTypes_P.ref == 'T') & (subTypes_P.alt == 'C')))
    mut_ACTG_P = sum(((subTypes_P.ref == 'A') & (subTypes_P.alt == 'C')) |
                    ((subTypes_P.ref == 'T') & (subTypes_P.alt == 'G')))
    sumMut_P = sum([mut_CTGA_P, mut_CGGC_P, mut_CAGT_P,
                mut_ATTA_P, mut_AGTC_P, mut_ACTG_P])


    # for recurrent without hypermutation
    mut_CTGA_R = sum(((subTypes_R.ref == 'C') & (subTypes_R.alt == 'T')) |
                    ((subTypes_R.ref == 'G') & (subTypes_R.alt == 'A')))
    mut_CGGC_R = sum(((subTypes_R.ref == 'C') & (subTypes_R.alt == 'G')) |
                    ((subTypes_R.ref == 'G') & (subTypes_R.alt == 'C')))
    mut_CAGT_R = sum(((subTypes_R.ref == 'C') & (subTypes_R.alt == 'A')) |
                    ((subTypes_R.ref == 'G') & (subTypes_R.alt == 'T')))
    mut_ATTA_R = sum(((subTypes_R.ref == 'A') & (subTypes_R.alt == 'T')) |
                    ((subTypes_R.ref == 'T') & (subTypes_R.alt == 'A')))
    mut_AGTC_R = sum(((subTypes_R.ref == 'A') & (subTypes_R.alt == 'G')) |
                    ((subTypes_R.ref == 'T') & (subTypes_R.alt == 'C')))
    mut_ACTG_R = sum(((subTypes_R.ref == 'A') & (subTypes_R.alt == 'C')) |
                    ((subTypes_R.ref == 'T') & (subTypes_R.alt == 'G')))
    sumMut_R = sum([mut_CTGA_R, mut_CGGC_R, mut_CAGT_R,
                mut_ATTA_R, mut_AGTC_R, mut_ACTG_R])

    # for recurrent with hypermutation
    mut_CTGA_HR = sum(((subTypes_HR.ref == 'C') & (subTypes_HR.alt == 'T')) |
                    ((subTypes_HR.ref == 'G') & (subTypes_HR.alt == 'A')))
    mut_CGGC_HR = sum(((subTypes_HR.ref == 'C') & (subTypes_HR.alt == 'G')) |
                    ((subTypes_HR.ref == 'G') & (subTypes_HR.alt == 'C')))
    mut_CAGT_HR = sum(((subTypes_HR.ref == 'C') & (subTypes_HR.alt == 'A')) |
                    ((subTypes_HR.ref == 'G') & (subTypes_HR.alt == 'T')))
    mut_ATTA_HR = sum(((subTypes_HR.ref == 'A') & (subTypes_HR.alt == 'T')) |
                    ((subTypes_HR.ref == 'T') & (subTypes_HR.alt == 'A')))
    mut_AGTC_HR = sum(((subTypes_HR.ref == 'A') & (subTypes_HR.alt == 'G')) |
                    ((subTypes_HR.ref == 'T') & (subTypes_HR.alt == 'C')))
    mut_ACTG_HR = sum(((subTypes_HR.ref == 'A') & (subTypes_HR.alt == 'C')) |
                    ((subTypes_HR.ref == 'T') & (subTypes_HR.alt == 'G')))
    sumMut_HR = sum([mut_CTGA_HR, mut_CGGC_HR, mut_CAGT_HR,
                    mut_ATTA_HR, mut_AGTC_HR, mut_ACTG_HR])

    tmp_list = []
    tmp_list.append(['Primary'] * 6 + ['Rec_noHM'] * 6 + ['Rec_HM'] * 6)
    tmp_frac = np.array(
        [mut_CTGA_P, mut_CGGC_P, mut_CAGT_P,
        mut_ATTA_P, mut_AGTC_P, mut_ACTG_P] +
        [mut_CTGA_R, mut_CGGC_R, mut_CAGT_R,
        mut_ATTA_R, mut_AGTC_R, mut_ACTG_R] +
        [mut_CTGA_HR, mut_CGGC_HR, mut_CAGT_HR,
        mut_ATTA_HR, mut_AGTC_HR, mut_ACTG_HR]
    )
    tmp_sum = np.array([sumMut_P] * 6 + [sumMut_R] * 6 + [sumMut_HR] * 6)
    tmp_list.append(tmp_frac / tmp_sum)
    tmp_list.append(['CTGA', 'CGGC', 'CAGT', 'ATTA', 'AGTC', 'ACTG'] * 3)
    df_HM_mutType = pd.DataFrame(np.array(tmp_list).T, columns=[
                                'sampleType', 'Fraction', 'MutType'])
    return(df_HM_mutType)

def cal_HM_SMratio(df_savi, savi_P, savi_R, savi_isHM):
    """
    """
    case = df_savi.CaseID.drop_duplicates().tolist()
    ratio_P = [np.NaN] * len(case)
    ratio_R = [np.NaN] * len(case)
    ratio_HR = [np.NaN] * len(case)

    for i in range(len(case)):  # for each patient
        tmp_P = df_savi[(df_savi.CaseID == case[i]) & (savi_P == 1)]
        tmp_R = df_savi[(df_savi.CaseID == case[i]) &
                        (savi_R == 1) & (savi_isHM == 0)]
        tmp_HR = df_savi[(df_savi.CaseID == case[i]) &
                        (savi_R == 1) & (savi_isHM == 1)]

        missense_P = sum(['MISSENSE' in strs for strs in tmp_P.Functional_Class])
        silent_P = sum(['SILENT' in strs for strs in tmp_P.Functional_Class])
        missense_R = sum(['MISSENSE' in strs for strs in tmp_R.Functional_Class])
        silent_R = sum(['SILENT' in strs for strs in tmp_R.Functional_Class])
        missense_HR = sum(['MISSENSE' in strs for strs in tmp_HR.Functional_Class])
        silent_HR = sum(['SILENT' in strs for strs in tmp_HR.Functional_Class])

        # filter out the samples that mutations fewer than 10
        if((len(tmp_P) > 10) & (missense_P > 0)):
            ratio_P[i] = silent_P / missense_P
        if((len(tmp_R) > 10) & (missense_R > 0)):
            ratio_R[i] = silent_R / missense_R
        if((len(tmp_HR) > 10) & (missense_HR > 0)):
            ratio_HR[i] = silent_HR / missense_HR

    tmp_list = []
    tmp_list.append(['Primary'] * len(ratio_P) + ['Rec_noHM'] * len(ratio_R)
                    + ['Rec_HM'] * len(ratio_HR))
    tmp_list.append(ratio_P + ratio_R + ratio_HR)
    tmp_list.append(['P'] * len(ratio_P) + ['R'] * len(ratio_R)
                    + ['R'] * len(ratio_HR))
    df_HM_ratio = pd.DataFrame(np.array(tmp_list).T, columns=[
                            'Groups', 'Ratio', 'PRmark'])
    return(df_HM_ratio)
