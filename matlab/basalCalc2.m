function v = basalCalc2(g)

global bas

v = (bas(6)-bas(1))/(bas(7)-bas(2))*(g-bas(2))*heaviside(g-bas(2))*heaviside(bas(7)-g)-...
    bas(1)/(bas(3)-bas(4))*(bas(3)-g)*heaviside((bas(3)-g))*heaviside(g-bas(4))+...
    bas(1)*heaviside(g-bas(4))+...
    (bas(6)-bas(1))*heaviside(g-bas(7));

end