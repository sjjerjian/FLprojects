function aMA = choose(p)

aMA = max(find([-eps cumsum(p)] < rand));

end