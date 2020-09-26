Unit Geral;

Interface
uses lu1;

Function Elevar(base : Float; expoente : Integer) : Float;

Function AbrirTexto(var arq : Text; nome : String) : Boolean;

Procedure ProcuraCar(s : String;var pont : Byte; car : Char; var code : Integer);
(*---------------------------------------------------------------------
  Proposito: Varre a cadeia s, comecando em pont, ate' encontrar o
             caracter car.
    Entrada: s - String[255]
             pont - posicao inicial da busca
             car - caracter procurado
      Saida: pont - posicao do caracter car
             code - codigo de retorno:
               0 - operacao concluida corretamente
               1 - pont inicial e' maior que o tamanho da cadeia
               2 - caracter nao encontrado
      Erros: Como previsto por code.
---------------------------------------------------------------------*)

Procedure ProcuraConjCar(s : String;var pont : Byte; car : String; var code : Integer);
(*---------------------------------------------------------------------
  Proposito: Varre a cadeia s, comecando em pont, ate' encontrar um dos
             caracteres de car.
    Entrada: s - String[255]
             pont - posicao inicial da busca
             car - caracteres procurados
      Saida: pont - posicao seguinte ao caracter encontrado
             code - codigo de retorno:
               0 - operacao concluida corretamente
               1 - pont inicial e' maior que o tamanho da cadeia
               2 - caracteres nao encontrados
      Erros: Como previsto por code.
---------------------------------------------------------------------*)

Procedure IgnoraCar(s : String;var pont : Byte; car : Char; var code : Integer);
(*---------------------------------------------------------------------
  Proposito: Varre a cadeia s, comecando em pont, enquanto encontrar o
             caracter car.
    Entrada: s - String[255]
             pont - posicao inicial da busca
             car - caracter a ser ignorado
      Saida: pont - posicao seguinte ao ultimo caracter car
             code - codigo de retorno:
               0 - operacao concluida corretamente
               1 - pont inicial e' maior que o tamanho da cadeia
               2 - ultimo caracter tambem e' car
      Erros: Como previsto por code.
---------------------------------------------------------------------*)

Procedure LeInt(s : String; var pont : Byte; var i, code : Integer);
(*---------------------------------------------------------------------
  Proposito: Varrer a cadeia s lendo um numero inteiro comecando em pont.
    Entrada: s - String[5]
             pont - posicao inicial da busca
      Saida: i - valor inteiro
             pont - posicao seguinte ao inteiro
             code - codigo de retorno:
               0 - Operacao concluida corretamente
               1 - Inteiro nao encontrado. Neste caso i vale 0
               2 - Inteiro excede faixa
               3 - Erro na conversao
      Erros: Como previsto por code
---------------------------------------------------------------------*)

Procedure LeReal(s : String; var pont: Byte; var j :Float ;var code :integer);
(*---------------------------------------------------------------------
  Proposito: Varrer a cadeia s lendo um numero real comecando em pont.
    Entrada: s - String[5]
             pont - posicao inicial da busca
      Saida: j - valor real
             pont - posicao seguinte ao inteiro
             code - codigo de retorno:
               0 - Operacao concluida corretamente
               1 - Point inicial maior que o tamanho da cadeia
               2 - Real nao encontrado. Neste caso j vale 0
               3 - Erro na conversao
      Erros: Como previsto por code
---------------------------------------------------------------------*)

Implementation

Function Elevar(base : Float; expoente : Integer) : Float;
var
  p : Float;
  i : Integer;
Begin
  p := 1.0;
  For i := 1 to expoente do p := p * base;
  Elevar := p;
End;

Function AbrirTexto(var arq : Text; nome : String) : Boolean;
var volta : Word;
Begin
  Assign(arq,nome);
  {$I-}
  Reset(arq);
  volta := IOResult;
  AbrirTexto := volta = 0;
  {$I+}
End;

Procedure ProcuraCar(s : String;var pont : Byte; car : Char; var code : Integer);
Begin
  code := 0;
  If pont <= Length(s) then
    begin
      While (s[pont] <> car) and (pont < Length(s)) do Inc(pont);
      If (pont = Length(s)) and (s[pont] <> car) then code := 2;
    end else code := 1;
End;

Procedure ProcuraConjCar(s : String;var pont : Byte; car : String; var code : Integer);
var
  encontrou : Boolean;
  posicao, i : Byte;
  c : String;
  menor_posicao : integer;
Begin
  s := s + ' ';
  code := 0;
  If pont <= Length(s) then
    begin
      c := Copy(s,pont,Length(s) - pont + 1);
      encontrou := False;
      menor_posicao := Length(c);  {menor_posicao indica o inicio do conjunto}
      For i := 1 to Length(car) do
        begin
          posicao := Pos(car[i],c);
          If (menor_posicao >= posicao) and (posicao > 0)
            then menor_posicao := posicao;
        end;
      If menor_posicao <> Length(c) then
        pont := pont + menor_posicao - 1 else code := 2;
    end else code := 1;
End;

Procedure IgnoraCar(s : String;var pont : Byte; car : Char; var code : Integer);
Begin
  code := 0;
  If pont <= Length(s) then
    begin
      While (s[pont] = car) and (pont < Length(s)) do Inc(pont);
      If (pont = Length(s)) and (s[pont] = car) then code := 2;
    end else code := 1;
End;

Procedure LeInt(s : String; var pont : Byte; var i, code : Integer);
var
  c : String;
  r : Float;
  codigo : Integer;
Begin
  code := 0;
  i := 0;
  c := '';
  If s[pont] in ['-','+'] then
    begin
      c := s[pont];
      Inc(pont);
    end;
  While (s[pont] in ['0'..'9']) and (pont <= Length(s)) do
    begin
      c := c + s[pont];
      Inc(pont);
    end;
  If c <> '' then
    begin
      Val(c,r,codigo);
      If codigo = 0 then
        If Abs(r) <= MAXINT then i := Trunc(r) else code := 2
        else code := 3;
    end else code := 1;
End;

Procedure LeReal(s : String; var pont: Byte; var j :Float ;var code :integer);
var aux_pont : Byte;
        code1: Integer;

Begin
  ProcuraConjCar(s,pont,'-+0123456789Ee.',code);
  aux_pont:=pont;
  ProcuraCar(s,aux_pont,' ',code1);
  If (code = 0) and (code1 = 0) then
    begin
      Val(Copy(s,pont,aux_pont - pont),j,code);
      If code <> 0 then code:=3;
    end;
  If (code = 0) and (code1=2) then
    Val(Copy(s,pont,Length(s)-pont+1),j,code);
End;

Begin
End.
