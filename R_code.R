#Mitter et al 2016

library(userfriendlyscience)
library(PMCMR)
library(Rcmdr)
library(coin)
options(scipen = 999)

#Fig_5
#Figure 5a
#WT
mq=c(476,455,589,455,509,613,878,788,574,470,522,543,501,619,644,529)
gus_naked=c(506.4507516,400.0501926,366.7057348,368.3613142,313.1467229,521.103247,545.8527425,473.4546329,414.2519787,447.1227202,414.0580917,398.8375488,380.5363008,379.0190864,419.3164789,397.337978)
gus_bioclay=c(434.5225697,537.3212576,572.1481697,404.0668319,489.4757722,537.6019621,542.9657486,543.7435215,431.0529028,442.6909771,423.6320503,361.6368202,406.6884351,379.4836076,364.0548899,352.4567172)
l4440=c(440.4312306,412.8825086,740.3003304,732.6993515,601.0695448,487.6049893,438.5686006,635.2338374,506.2871921,485.4731696,360.7145246,408.8430045,447.1101129,460.805071,622.7308165,657.8322937)
two_b=c(437.3078278,428.9133046,483.2695802,548.2762745,424.5072598,655.8430879,673.8530287,780.793623,413.2287663,557.6841135,454.8120719,628.6126026,728.8983555,431.2129426,469.4434101,466.8503749)
ldh=c(865.5848031,546.3712308,558.0085599,442.8708797,715.2306024,531.1353423,369.3582536,617.7371673,661.2298057,474.964126,468.5175731,479.8913825,462.449421,626.7318906,640.406223,498.2266169)

all_activity = c(mq, gus_naked, gus_bioclay, l4440, two_b, ldh)

treatment=c(rep("mq",
                length(mq)),
            rep("gus_naked",
                length(gus_naked)),
            rep("gus_bioclay",
                length(gus_bioclay)),
            rep("l4440",
                length(l4440)),
            rep("two_b",
                length(two_b)),
            rep("ldh",
                length(ldh))
)
experiment=data.frame(all_activity,treatment)

#test for normality (threshold = 0.01 given the tests are quite robust)
shapiro.test(mq)
shapiro.test(gus_naked)
shapiro.test(gus_bioclay)
shapiro.test(l4440)
shapiro.test(two_b)
shapiro.test(ldh)

#null hypothesis -> all groups normal - NOT REJECTED

#test for equal variances
leveneTest(all_activity~treatment, data=experiment)
#null hypothesis -> variances same for all groups - NOT REJECTED

#ANOVA then pairwise t test with holm adjustment
fit <- aov(all_activity~treatment, data=experiment)
pairwise.t.test(all_activity, treatment, p.adj = "holm")

#Figure 5b
#rdr6
mq=c(666.5920573,798.7708326,534.8607055,438.0662022,778.156637,434.0394856,511.920889,409.0039689,703.9472216,849.6336698,493.794251,340.9954782,476.0063279,586.6278142,578.0933978,582.2528012)
gus_naked=c(685.1519238,586.2033966,479.7860186,642.913244,597.9235745,503.480292,464.897759,579.500643,546.5996762,499.0580109,482.825204,388.4233049,535.3115644,505.4821381,641.5349558,524.7302206)
gus_bioclay=c(508.2101485,700.5133652,447.4379401,605.8279377,333.4488746,444.4496665,652.6108734,700.0683604,533.0747524,529.6836837,450.122061,524.4442716,488.7597367,497.3621127,490.3120412,633.6793683)
l4440=c(571.3947454,536.9072912,532.2696333,413.3470774,467.2919079,397.0003871,571.2685329,805.769811,479.4042203,442.3843028,493.5510955,771.4799715,414.2711082,470.4928394,488.5101947,485.7356925)
two_b=c(556.8044942,546.5144189,461.4792126,688.2985497,487.600733,542.1384929,610.7052034,689.4803942,456.3213336,434.6173841,489.2068345,435.5745154,519.4771682,536.5032039,511.6647207,445.7231717)
ldh=c(467.647423,531.4474401,723.8934533,692.472472,550.5845424,618.6863964,613.4382178,523.9523857,645.6217627,524.747695,423.69807,427.7954483,448.9750848,420.1292172,699.2412305,527.8256575)

all_activity = c(mq, gus_naked, gus_bioclay, l4440, two_b, ldh)

treatment=c(rep("mq",
                length(mq)),
            rep("gus_naked",
                length(gus_naked)),
            rep("gus_bioclay",
                length(gus_bioclay)),
            rep("l4440",
                length(l4440)),
            rep("two_b",
                length(two_b)),
            rep("ldh",
                length(ldh))
)
experiment=data.frame(all_activity,treatment)

#test for normality (threshold = 0.01 given the tests are quite robust)
shapiro.test(mq)
shapiro.test(gus_naked)
shapiro.test(gus_bioclay)
shapiro.test(l4440)
shapiro.test(two_b)
shapiro.test(ldh)

#null hypothesis -> all groups normal - Rejected for L4440, but not essential for pairwise comparison

#test for equal variances
leveneTest(all_activity~treatment, data=experiment)
#null hypothesis -> variances same for all groups - NOT REJECTED

#One way test with equal variances - holm multiple correction
oneway(x=experiment$treatment,y=experiment$all_activity, levene=FALSE, posthoc="holm")

#ANOVA then holm adjustment
fit <- aov(all_activity~treatment, data=experiment)
pairwise.t.test(all_activity, treatment, p.adj = "holm")

#Fig_6
#Figure 6a

#1 day results

#Vectors for each treatment
ldh_lesions=c(74,14,92,68,74,11,75,34)
dsRNA_lesions=c(7,11,27,7,3,15,5,2,14,18,12,4,2,7,4,5)
bioclay_lesions=c(14,13,4,27,9,19,6,4,7,5,5,13,7,15,8,5)

# ldh_lesions=c(148,25,167,102)
# dsRNA_lesions=c(21,29,39,11,5,22,9,7)
# bioclay_lesions=c(21,18,9,40,16,34,14,9)

#Convert to a data frame
lesions=c(ldh_lesions,dsRNA_lesions,bioclay_lesions)
treatment=c(rep("ldh",
                length(ldh_lesions)),
            rep("dsrna",
                length(dsRNA_lesions)),
            rep("bioclay",
                length(bioclay_lesions)))
experiment=data.frame(lesions,treatment)


#test for normality (threshold = 0.01 given the tests are quite robust)
shapiro.test(ldh_lesions)
shapiro.test(dsRNA_lesions)
shapiro.test(bioclay_lesions)
#null hypothesis -> all groups normal - NOT REJECTED

#test for equal variances
leveneTest(lesions~treatment, data=experiment)
#null hypothesis -> variances same for all groups - REJECTED

#Test for difference in means
kruskal.test(lesions~treatment, data=experiment)
#null hpypothesis --> means of all groups the same --> REJECTED

#Post hoc test
posthoc.kruskal.nemenyi.test(x=experiment$lesions, g=experiment$treatment, dist="Chisquare")


#5 day results

ldh_lesions = c(57,78,113,65,112,33,82,67,57,55,68,56,120,43,73,90)
dsRNA_lesions = c(54,23,39,32,31,14,29,12,93,23,39,24,31,14,18,8)
bioclay_lesions = c(35,13,31,51,17,33,8,31,63,19,16,47,33,17,18,18)

# ldh_lesions = c(114,133,181,121,232,76,155,157)
# dsRNA_lesions = c(147,46,78,56,62,28,47,20)
# bioclay_lesions = c(98,32,47,98,50,50,26,49)


#Convert to a data frame
lesions=c(ldh_lesions,dsRNA_lesions,bioclay_lesions)
treatment=c(rep("ldh",
                length(ldh_lesions)),
            rep("dsrna",
                length(dsRNA_lesions)),
            rep("bioclay",
                length(bioclay_lesions)))
experiment=data.frame(lesions,treatment)

#boxplot
plot(lesions~treatment, data=experiment)

#test for normality (threshold = 0.01 given the tests are quite robust)
shapiro.test(ldh_lesions)
shapiro.test(dsRNA_lesions)
shapiro.test(bioclay_lesions)
##null hypothesis -> all groups normal - REJECTED, dsRNA p =0.003

#test for equal variances
bartlett.test(lesions~treatment, data=experiment)
leveneTest(lesions~treatment, data=experiment)
#null hypothesis --> all groups have same variance NOT REJECTED

#Therefore use Kruskal-Wallis non-parametric test, as distributions are all right-skewed (i.e. similar)

#Test for difference in means
kruskal.test(lesions~treatment, data=experiment)
#null hpypothesis --> means of all groups the same --> REJECTED

#Post hoc test - non parametric using chisquare distibution due to tied data
posthoc.kruskal.nemenyi.test(x=experiment$lesions, g=experiment$treatment, dist="Chisquare")

#Figure 6b

#Day 5 
water_lesions = c(12,16,28,41,18,81,33,30,37,88,63,30,22,23,33,13,29,39,42)
ldh_lesions = c(38,40,48,41,53,15,11,15,17,25,21,7,25,31,19)
dsRNA_lesions = c(3,3,1,1,3,1,0,4,1,2,6,0,7,2,0,2,5,1,4,2)
bioclay_lesions = c(1,4,1,3,1,1,2,0,1,3,1,9,1,1,2,6,6,1,1,0,2,0,2,7,1)

#convert to data frame
lesions=c(water_lesions,ldh_lesions,dsRNA_lesions,bioclay_lesions)
treatment=c(rep("water",
                length(water_lesions)),
            rep("ldh",
                length(ldh_lesions)),
            rep("dsrna",
                length(dsRNA_lesions)),
            rep("bioclay",
                length(bioclay_lesions)))

experiment=data.frame(lesions,treatment)


#test for normality (threshold = 0.01 given the tests are quite robust)
shapiro.test(water_lesions)
shapiro.test(ldh_lesions)
shapiro.test(dsRNA_lesions)
shapiro.test(bioclay_lesions)
#null hypothesis -> all groups normal - REJECTED

#test for equal variances
bartlett.test(lesions~treatment, data=experiment)
levene.test(lesions~treatment, data=experiment)
#null hypothesis -> variances same for all groups - REJECTED


#Test for difference in population distributions using Kruskal-Wallis due to non-normal distribution
kruskal.test(lesions~treatment, data=experiment)
#null hpypothesis --> medians of all groups the same --> REJECTED

#Post hoc test - non parametric using chisquare distibution due to tied data
posthoc.kruskal.nemenyi.test(x=experiment$lesions, g=experiment$treatment, dist="Chisquare")

#Day 20
water_lesions = c(31,30,29,25,17,13,15,17,16,31)
ldh_lesions = c(31,18,33,3,12,11,14,13,3,16,21,2,15)
dsRNA_lesions = c(8,20,17,41,23,24,19,28,29,26)
bioclay_lesions = c(0,10,1,2,3,1,3,4,5,7,3,1,2,1,4)

#convert to data frame
lesions=c(water_lesions,ldh_lesions,dsRNA_lesions,bioclay_lesions)
treatment=c(rep("water",
                length(water_lesions)),
            rep("ldh",
                length(ldh_lesions)),
            rep("dsrna",
                length(dsRNA_lesions)),
            rep("bioclay",
                length(bioclay_lesions)))
experiment=data.frame(lesions,treatment)

#boxplot
plot(lesions~treatment, data=experiment)

#test for normality (threshold = 0.01 given the tests are quite robust)
shapiro.test(water_lesions)
shapiro.test(ldh_lesions)
shapiro.test(dsRNA_lesions)
shapiro.test(bioclay_lesions)
#null hypothesis -> all groups normal - NOT REJECTED, though some <0.05

#test for equal variances
bartlett.test(lesions~treatment, data=experiment)
leveneTest(lesions~treatment, data=experiment)
#null hypothesis -> variances same for all groups - REJECTED

#Test for difference in population distributions using Kruskal-Wallis due to non-normal distribution
kruskal.test(lesions~treatment, data=experiment)
#null hpypothesis --> medians of all groups the same --> REJECTED

#Post hoc test - non parametric using chisquare distibution due to tied data
posthoc.kruskal.nemenyi.test(x=experiment$lesions, g=experiment$treatment, dist="Chisquare")

#Fig_7
#Figure 7a

ldh=c(100,62.5,70,66.7)
dsrna=c(56.25,50,90,86.11)
bioclay=c(6.25,25,30,36.11)

all_percent = c(ldh, dsrna, bioclay)

treatment=c(rep("ldh",
                length(ldh)),
            rep("dsrna",
                length(dsrna)),
            rep("bioclay",
                length(bioclay))
)

experiment=data.frame(all_percent,treatment)

#boxplot
plot(all_percent~treatment, data=experiment)

#test for normality (threshold = 0.01 given the tests are quite robust)
shapiro.test(ldh)
shapiro.test(dsrna)
shapiro.test(bioclay)

#null hypothesis -> all groups normal - NOT REJECTED

#test for equal variances
leveneTest(all_percent~treatment, data=experiment)
#null hypothesis -> variances same for all groups - NOT REJECTED

#One way test with equal variances - holm multiple correction
oneway(x=experiment$treatment,y=experiment$all_percent, levene=FALSE, posthoc="holm")

#ANOVA then holm adjustment
fit <- aov(all_percent~treatment, data=experiment)
pairwise.t.test(all_percent, treatment, p.adj = "holm")

#Sup_Fig_3

mq=c(509,530,490,582)
depc=c(374,453,393,380)
ldh_300=c(427,518,449,572)
ldh_900=c(543,459,428,412)
L4440=c(394,519,415,390)
two_b=c(522,442,417,505)
col=c(5,11,10)
gus_100=c(490,367,469,448)
gus_200=c(442,444,309,458)
gus_300=c(446,307,427,452)
bioclay_100=c(458,425,295,347)
bioclay_200=c(356,390,382,399)
bioclay_300=c(391,274,313,363)

all_gus=c(mq,depc,ldh_300,ldh_900,L4440,two_b,col,
          gus_100,gus_200,gus_300,bioclay_100,
          bioclay_200,bioclay_300)

treatment=c(rep("mq",
                length(mq)),
            rep("depc",
                length(depc)),
            rep("ldh_300",
                length(ldh_300)),
            rep("ldh_900",
                length(ldh_900)),
            rep("L4440",
                length(L4440)),
            rep("two_b",
                length(two_b)),
            rep("col",
                length(col)),
            rep("gus_100",
                length(gus_100)),
            rep("gus_200",
                length(gus_200)),
            rep("gus_300",
                length(gus_300)),
            rep("bioclay_100",
                length(bioclay_100)),
            rep("bioclay_200",
                length(bioclay_200)),
            rep("bioclay_300",
                length(bioclay_300))
)
experiment=data.frame(all_gus,treatment)

#test for normality (threshold = 0.01 given the tests are quite robust)
shapiro.test(mq)
shapiro.test(depc)
shapiro.test(ldh_300)
shapiro.test(L4440)
shapiro.test(two_b)
shapiro.test(col)
shapiro.test(gus_100)
shapiro.test(gus_200)
shapiro.test(gus_300)
shapiro.test(bioclay_100)
shapiro.test(bioclay_200)
shapiro.test(bioclay_300)
##null hypothesis -> all groups normal - NOT REJECTED, though <0.05

#test for equal variances
leveneTest(all_gus~treatment, data=experiment)
#null hypothesis -> variances same for all groups - NOT REJECTED

#ANOVA then holm adjustment
fit <- aov(all_gus~treatment, data=experiment)
summary(fit)
pairwise.t.test(all_gus, treatment, p.adj = "holm")

#Sup_Fig_5

#TRIAL 3

water_only=c(263,342,317,208,591,470,267,99,100,172,226,367,254)
ir54=c(307,208,314,127,163,181,218,90,348,57,60,91,41)
dsRNA_2b=c(106,55,23,38,49,45,75,28,40,11,31,23,32)



all_treatments=c(water_only,ir54,dsRNA_2b)

treatment=c(rep("water_only",
                length(water_only)),
            rep("ir54",
                length(ir54)),
            rep("dsRNA_2b",
                length(dsRNA_2b))
            )

experiment=data.frame(all_treatments,treatment)


#test for normality (threshold = 0.01 given the tests are quite robust)
shapiro.test(water_only)
shapiro.test(ir54)
shapiro.test(dsRNA_2b)

#test for equal variances
#bartlett.test(all_gus~treatment, data=experiment)
leveneTest(all_treatments~treatment, data=experiment)
#null hypothesis -> variances same for all groups - NOT REJECTED

#Due to unequal variances:
kruskal.test(all_treatments~treatment, data=experiment)
posthoc.kruskal.nemenyi.test(x=experiment$all_treatments, g=experiment$treatment, dist="Chisquare")


#TRIAL 4

water_only=c(11,45,32,48,56,45,61,9)
ir54=c(30,14,52,27,60,48,41,37)
dsRNA_2b=c(5,15,9,10,8,11,10,7)
GUS=c(23,43,36,47,32,16,34,42)

all_treatments=c(water_only,ir54,dsRNA_2b, GUS)

treatment=c(rep("water_only",
                length(water_only)),
            rep("ir54",
                length(ir54)),
            rep("dsRNA_2b",
                length(dsRNA_2b)),
            rep("GUS",
                length(GUS))
)

experiment=data.frame(all_treatments,treatment)

#boxplot
plot(all_treatments~treatment, data=experiment)

#test for normality (threshold = 0.01 given the tests are quite robust)
shapiro.test(water_only)
shapiro.test(ir54)
shapiro.test(dsRNA_2b)
shapiro.test(GUS)

#test for equal variances
#bartlett.test(all_gus~treatment, data=experiment)
leveneTest(all_treatments~treatment, data=experiment)
#null hypothesis -> variances same for all groups - NOT REJECTED

#Due to unequal variances:
kruskal.test(all_treatments~treatment, data=experiment)
posthoc.kruskal.nemenyi.test(x=experiment$all_treatments, g=experiment$treatment, dist="Chisquare")
