from calculateMetrics import *

class DiffCI(BasePara):
    #compare two metrices need normalization
    def __init__(self,path,control_path,resolution,chromosome,out_name="DRF",
                diffCI_size=300000):
        super().__init__(path,resolution,chromosome,out_name)
        self.path = path
        self.control_path = control_path
        treat = loadWithNorm(path,log = True).values
        control = loadWithNorm(control_path,log = True).values
        smooth = 2
        self.treat = ndimage.median_filter(treat,smooth)
        self.control = ndimage.median_filter(control,smooth)
        self.matrix = self.treat-self.control
        self.diffCI_size = diffCI_size

        if self.treat.shape[0] != self.control.shape[0]:
            print("Error: Input/control matrix have different shape")
            exit(1)

    def getDiffCI(self):
        CI_Bin = round(self.diffCI_size/self.resolution)

        array = np.zeros(self.matrix_shape)
        for i in range(self.matrix_shape):
            if(i - CI_Bin < 0 or i + CI_Bin >= self.matrix_shape): continue
            matA = self.matrix[i-CI_Bin:i,i-CI_Bin:i]
            matB = self.matrix[i+1:i+CI_Bin+1,i+1:i+CI_Bin+1]
            A = np.triu(matA,1).sum()
            B = np.triu(matB,1).sum()
            C = self.matrix[i-CI_Bin: i, i+1: i+CI_Bin+1].sum()
            if np.isnan(A) or np.isnan(B) or np.isnan(C): continue  #skip NaN value
            S = 2*C-(A+B)
            if S < 0:
                array[i] = -np.log1p(-S)
            else:
                array[i] = np.log1p(S)

        return super().makeDF(array,"diffCI")

    #def getCSV(self):
    #    super().makeCSV(self.getDRF())
