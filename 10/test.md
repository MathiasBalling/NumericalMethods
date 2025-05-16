# Test case for Trapezoidal

Pga. jeres problemer med at implementere Trapez-metoden i Mandatory 4 får I her et testproblem som er hugget fra kemisk kinetik:

y1'(t)= - 0.05 y1(t) + 10000 y2(t)y3(t); y1(0)=1
y2'(t)= 0.05 y1(t) - 10000 y2(t)y3(t) - 30000000y2(t)^2; y2(0)=0
y3'(t)=30000000y2(t)^2; y3(0)=0

Brug Trapez-metoden med h=0.005 og find y1(5),y2(5),y3(5) (altså svarende til 1000 skridt). I bør få ca. {0.8713, 0.000022, 0.1286} (medmindre jeg har lavet en fejl hvilket man aldrig kan udelukke :)). Hvis I ikke har implementeret Trapezmetoden korrekt bliver dette meget tydeligt når I forsøger at løse dette problem.
Skriv kommentar
