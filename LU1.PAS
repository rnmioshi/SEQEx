Unit LU1;

Interface

const
  TNNearlyZero = 1E-35;
  TNArraySize = 30;               { Size of the matrix }

type
  Float = Extended; { 8 byte real, requires 8087 math chip }
  TNvector = array[1..TNArraySize] of Float;
  TNmatrix = array[1..TNArraySize] of TNvector;
  glindx   = Array[1..TNArraySize] of Integer;


Procedure Ludcmp ( var a: TNmatrix; n: integer;
                   var indx: glindx; var d: Float;var erro : Byte);

Procedure Lubksb (a: TNmatrix; n: integer; indx: glindx;b: TNvector; var d: TNvector);

Implementation

Procedure Ludcmp ( var a: TNmatrix; n: integer;
                   var indx: glindx; var d: Float;var erro:Byte);

Const
   tiny=1.0e-20;
var
  k,j,imax,i : integer;
  sum,dum,big: Float;
           vv: TNvector;

BEGIN
  erro:=0;
  d:=1.0;
  For i:= 1 to n do
    begin
      big:=0.0;
      for j:= 1 to n do if (abs(a[i,j])>big) then big:= abs(a[i,j]);
      if (big = 0.0) then
      begin
        writeln('Pause in LUDCMP - singular matrix');
        erro:=1;
      end;
      vv[i]:=1.0/(big);
    end;

  For j:= 1 to n do
    begin

      For i:= 1 to (j-1) do
        begin
           sum:=a[i,j];
           For k:=1 to i-1 do
            begin
              sum:= sum-a[i,k]*a[k,j];
            end;
          a[i,j]:= sum
        end;

      big:= 0.0;

      For i:= j to n do
        begin
          sum:= a[i,j];
          For k:= 1 to j-1 do
            begin
              sum:= sum-a[i,k]*a[k,j];
            end;
          a[i,j]:= sum;
          dum:= vv[i] * abs(sum);
          If (dum >= big) then
            begin
              big:= dum;
              imax:= i;
            end;
        end;
      If (j <> imax) then
        begin
          For k:= 1 to n do
            begin
              dum:= a[imax,k];
              a[imax,k]:= a[j,k];
              a[j,k]:= dum
            end;
          d:=-d;
          vv[imax]:=vv[j]
        end;
      indx[j]:= imax;
      If (a[j,j] = 0.0) then a[j,j]:= tiny;
      If (j <> n) then
        begin
          dum:= 1.0/a[j,j];
          For i:= j+1 to n do
          begin
            a[i,j]:= a[i,j]*dum;
          end
        end;
    end;
END;

Procedure Lubksb(a: TNmatrix; n: integer; indx: glindx;b: TNVector;var d: TNVector);
var j,ip,l,i : integer;
    sum:       Float;
BEGIN
  l:=0;
  For i:= 1 to n do
  begin
    ip:= indx[i];
    sum:= b[ip];
    b[ip]:= b[i];
    If ((l) <> 0) then
    begin
      For j:= l to i-1 do
      begin
        sum:= sum-a[i,j]*b[j];
      end;
    end
    else If (sum <> 0.0) then
    begin
      l:= i;
    end;
    b[i]:= sum;
  end;
  For i:= n DOWNTO 1 do
  begin
    sum:= b[i];
    If (i<n) then
    begin
        For j:= i+1 to n do
        begin
          sum:= sum-a[i,j]*b[j];
        end;
    end;
    b[i]:= sum/a[i,i]
  end;
  d:=b;
END;

END.