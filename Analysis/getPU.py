fileName = "/users/cl7359/trigger_studies/CMSSW_5_0_1/src/RecoLuminosity/LumiDB/getLumi_out_179828_v5_Corr.txt"
f1 = open(fileName, 'r')

temp = f1.read()
temp = temp.split("\n")

for i in range(len(temp)):
	if ((i%3) == 0):
		print temp[i], temp[i+2]