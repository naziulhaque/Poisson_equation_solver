function x=conc(a,b)

for i=1:length(a)
    x(i)=a(i);
end
for i=1:length(b)
    x(i+length(a))=b(i);
end

end
