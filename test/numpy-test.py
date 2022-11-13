import numpy as np

# weird cases that produce a root that doesn't evaluate to zero.
#
#  case 90
# (-140.66232642249253+0j)  --> 3.8840212075494556
# (-0.8495747907453617-0.18984925426537902j)  --> 4.895797037288796e-13
# (-0.8495747907453617+0.18984925426537902j)  --> 4.895797037288796e-13
# (0.014327558420734632-0.9296444676629798j)  --> 5.176163245255504e-13
# (0.014327558420734632+0.9296444676629798j)  --> 5.176163245255504e-13
# (0.055116133964008746+0j)  --> 0.0
# (0.6838441575682882+0j)  --> 1.1048939541069558e-12

#  case 113
# (-0.8104330345121576+0j)  --> 9.094947017729282e-13
# (-0.7930272777739451-0.7498282775774507j)  --> 4.004940581991954e-12
# (-0.7930272777739451+0.7498282775774507j)  --> 4.004940581991954e-12
# (0.020968472437567498-1.0338747234876766j)  --> 4.466550987380886e-12
# (0.020968472437567498+1.0338747234876766j)  --> 4.466550987380886e-12
# (0.6475070012645461-0.6270937740733488j)  --> 2.332658883538351e-12
# (0.6475070012645461+0.6270937740733488j)  --> 2.332658883538351e-12
# (1.9014675801855059+0j)  --> 3.431068762438372e-10
# (28.257072684212886+0j)  --> 0.2127890375388688

print('\n case 90')

c = np.flip(np.array([   6.60460233688354,
                         935.171142578125,
                         867.901550292969,
                         352.381774902344,
                         320.264373779297,
                        -332.592407226562,
                        -398.892456054688,
                         22.9384136199951 ]))

roots = np.polynomial.polynomial.polyroots(c)

for r in roots:
    val = np.abs(np.polynomial.polynomial.polyval(r, c))
    print(r, ' -->', val)


print('\n case 113')
c = np.flip(np.array([  19.7673424520622,
                        -575.209969604783,
                        454.342995592928,
                        422.066614439551,
                        781.088236436477,
                        358.472540727891,
                        788.011551825783,
                        14.8576676207117,
                        330.581477423860,
                        890.813059508382     ]))

roots = np.polynomial.polynomial.polyroots(c)

for r in roots:
    val = np.abs(np.polynomial.polynomial.polyval(r, c))
    print(r, ' -->', val)
