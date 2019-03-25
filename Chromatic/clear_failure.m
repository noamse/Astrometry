function astcat = clear_failure(astcat)
Flag = true(size(astcat));
for i=1:numel(astcat)
    if (numel(fields(astcat(i).UserData.R))<=2)
        Flag(i)=false;
    end
end

astcat = astcat(Flag);




end