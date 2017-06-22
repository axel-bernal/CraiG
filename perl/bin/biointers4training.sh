#locs1 locs2 ll|ii|cc                                                          
biointers.py $1 $2 -by $3 > a$$;
biodiff.py $1 a$$ -by $3 > b$$;
biodiff.py $2 a$$ -by $3 > c$$;
biodiff.py $1 b$$ -by cc | biodiff.py - c$$ -by cc;
rm a$$; rm b$$; rm c$$;
