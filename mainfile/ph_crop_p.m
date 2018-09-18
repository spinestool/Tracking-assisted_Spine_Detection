function [ ph_crop_pozitive]=ph_crop_p(ph_crop)

[m1,n1]=size(ph_crop);ph_crop_pozitive=zeros(m1,n1);
           for ph_in1=1:m1
            for ph_in2=1:n1
                if ph_crop(ph_in1,ph_in2)<=0;
                    ph_crop_pozitive(ph_in1,ph_in2)=1;
                end
            end
           end
           