program MauchSpg;
{ Maucha R (1932) Hydrochemische Metoden in der Limnologie. Binnengewasser
    12, 173p.
  Broch E S & Yake W (1969) A modification of Maucha's ionic diagram
    to include ionic concentrations. Limnology & Oceanography 14 p 933-935 }
{ M J Silberbauer, Zoology University of Cape Town, September 1989 }
{ Modified to write to an Arc/Info Ungenerate file - MJS October 1992 }

{$R+} { check ranges }
Uses Graph;
var
   ChemFile   : Text;
   SpgFile    : Text;
   ChemValid  : boolean;
   Acidity    : boolean;
   SpgName    : string [64];
   VleiCode   : string [10];
   TextItem   : string [20];
   Pad        : string [1];
   X          : real;
   Y          : real;
   XAspect    : word;
   YAspect    : word;
   Aspect     : real;
   TotalArea  : real;
   IonArea    : real;
   TotalRadius: real;
   IonRadius  : real;
   XCentre    : integer;
   YCentre    : integer;
   SeedX      : integer;
   SeedY      : integer;
   XScreen    : integer;
   YScreen    : integer;
   I          : integer;
   Angle      : integer;
   Scale      : integer;
   StartAngle : real;
   IonAngle   : real;
   EndAngle   : real;
   Equivalents: real;
   EquivSum   : real;
   Vertex     : array[1..5,1..2] of integer;
   ChemInput  : array[1..9] of real;
   Gd, Gm     : integer;
   Temperature, Salinity, Conductivity, pH , TDS, TSS : real;
   CO3        , HCO3    , H2CO3       , Ca , Mg , Na  : real;
   K          , Cl      , S           , SO4, PO4      : real;
   SiO4       , TP      , KN          , NO3, NO2, NH4 : real;
   Polyphenol , WaterColour                           : real;
   Latitude   : real;
   Longitude  : real;

function r (Degree : real) : real;
begin
   r := Degree * Pi / 180;
end;

function ScaleX (XValue : real) : integer;      { scales x values }
begin
   XValue := XValue + (GetMaxX div 2);
   ScaleX := Round (XValue);
end;

function ScaleY (YValue : real) : integer;      { scales y values }
begin
   YValue := YValue * Aspect;
   YValue := YValue + (GetMaxY div 2);
   ScaleY := Round (YValue);
end;

procedure Equilivate;            { Converts mg/litre to meq/litre }
begin
   ChemInput[06] := CO3*2 / 50;  { "mg/l" originally obtained by multiplying}
   ChemInput[05] := HCO3  / 50;  {  meq/l by 50 (« mol. wt. of CaCO3)       }
   ChemInput[09] := H2CO3 / 50;
   ChemInput[01] := Ca*2  / 40.08;
   ChemInput[02] := Mg*2  / 24.31;
   ChemInput[08] := Na    / 22.99;
   ChemInput[07] := K     / 39.10;
   ChemInput[04] := Cl    / 35.45;
   ChemInput[03] := S*2   / 32.07;
end;


begin
   Assign (ChemFile, 'c:\data\vlei\CHEMDATA.DAT');
   Reset (ChemFile);
   Gd := Detect;
   InitGraph (Gd, Gm, '');
   if GraphResult <> grOK then Halt(1);
   SetTextStyle (TriplexFont, HorizDir, 1);
   SetTextJustify (LeftText, TopText);
   GetAspectRatio (XAspect, YAspect);
   Aspect := (XAspect / YAspect) * 0.9;

   Repeat
      For I := 1 to 9 do
      begin
         ChemInput[I] := -1;
      end;
      ReadLn (ChemFile, Latitude, Longitude, Latitude, Longitude,
              Temperature, Salinity, Conductivity, pH, TDS, TSS,
              CO3, HCO3, H2CO3, PO4, SiO4, TP, KN, NO3, NO2, NH4,
              Ca, Mg, Na, K, Cl, S,
              Polyphenol, WaterColour, Pad, Pad, VleiCode);
      Equilivate;

      ChemValid := true;
      EquivSum := 0;
      For I := 1 to 9 do
      begin
         If (ChemInput[I] < 0)  then ChemValid := false;
         EquivSum := EquivSum + ChemInput[I];
      end;

      Acidity := false;           { check for special case of acidity }
      if (H2CO3 > 0.01) then
      begin
         Acidity := true;
         ChemInput[05] := ChemInput[09]; { make bicarb represent acid.. }
      end;

      If (ChemValid) and (EquivSum > 0) then
      begin
         SpgName := Concat ('C:\DATA\VLEI\', VleiCode);
         WriteLn ('Opening ', SpgName);
         Assign (SpgFile, SpgName);
         Rewrite (SpgFile);

         XCentre     := 0;
         YCentre     := 0;
         TotalArea   := 4000 * Ln (EquivSum+1);
         TotalRadius := Sqrt (0.125 * TotalArea / Sin(r(22.5)));
         SetFillStyle (EmptyFill, 0);
         Ellipse     (ScaleX (XCentre), ScaleY (YCentre), 0, 360,
                      Round(TotalRadius), Round(TotalRadius * Aspect));

         FillEllipse (ScaleX (XCentre), ScaleY (YCentre),
                      Round(TotalRadius), Round(TotalRadius * Aspect));

         for Angle := 0 to 7 do
         begin
            Equivalents := ChemInput [Angle + 1];
            IonArea     := (Equivalents / EquivSum) * TotalArea;
            IonRadius   := IonArea / (TotalRadius * Sin(r(22.5)));
            StartAngle  := Angle * 45;
            IonAngle    := StartAngle + 22.5;
            EndAngle    := (Angle + 1) * 45;
            Vertex[1,1] := ScaleX (XCentre);
            Vertex[1,2] := ScaleY (YCentre);
            Vertex[2,1] := ScaleX (TotalRadius * Cos(r(StartAngle)));
            Vertex[2,2] := ScaleY (TotalRadius * Sin(r(StartAngle)));
            Vertex[3,1] := ScaleX (IonRadius   * Cos(r(IonAngle  )));
            Vertex[3,2] := ScaleY (IonRadius   * Sin(r(IonAngle  )));
            Vertex[4,1] := ScaleX (TotalRadius * Cos(r(EndAngle  )));
            Vertex[4,2] := ScaleY (TotalRadius * Sin(r(EndAngle  )));
            Vertex[5,1] := ScaleX (XCentre);
            Vertex[5,2] := ScaleY (YCentre);
            SeedX       := ScaleX (IonRadius/2.0 * Cos(r(IonAngle  )));
            SeedY       := ScaleY (IonRadius/2.0 * Sin(r(IonAngle  )));
            SetFillStyle (SolidFill, Angle+2);
            FillPoly (5, Vertex);
            SetColor (0);
            SetLineStyle (SolidLn, 0, ThickWidth);
                                      { Write to the Arc/Info Ungenerate file }
            WriteLn (SpgFile, Angle + 1, ',', SeedX, ',', SeedY);
            WriteLn (SpgFile, Vertex[1,1], ',', Vertex[1,2]);
            WriteLn (SpgFile, Vertex[2,1], ',', Vertex[2,2]);
            WriteLn (SpgFile, Vertex[3,1], ',', Vertex[3,2]);
            WriteLn (SpgFile, Vertex[4,1], ',', Vertex[4,2]);
            WriteLn (SpgFile, Vertex[5,1], ',', Vertex[5,2]);
            WriteLn (SpgFile, 'END');
            Line
            (Vertex[1,1], Vertex[1,2], Vertex[2,1], Vertex[2,2]);
            Line
            (Vertex[4,1], Vertex[4,2], Vertex[5,1], Vertex[5,2]);
            SetLineStyle (SolidLn, 0, NormWidth);
            SetColor (Angle+2);
         end;
         WriteLn (SpgFile, 'END');
         Close (SpgFile);
         X := TotalRadius * Sin(r(292.5));
         Y := TotalRadius * Cos(r(292.5));
         SetUserCharSize (2, 4, 2, 5);
         SetTextJustify (RightText, BottomText);
         Str (Round(TDS), TextItem);
         OutTextXY (GetMaxX-10, GetMaxY-10, 'TDS = ' + TextItem + ' mg/l');
         SetTextJustify (LeftText, BottomText);
         If (Acidity) then OutTextXY (20, GetMaxY-10, 'Acid');
         SetUserCharSize (3, 1, 3, 2);
         SetTextJustify (LeftText, TopText);
         OutTextXY (10,10,VleiCode);
         ReadLn;
         ClearDevice;
      end;
   until Eof (ChemFile);
   Close (ChemFile);
   SetUserCharSize (4, 3, 3, 3);
   EquivSum := 1;
   XCentre := 150;
   YCentre := 0;
   OutTextXY (XCentre + 40, YCentre, '[Equivalents]');
   SetTextJustify (LeftText, CenterText);
   for I := 1 to 5 do
   begin
      YCentre     := YCentre + 50;
      TotalArea   := 4000 * Ln (EquivSum+1);
      TotalRadius := Sqrt (0.125 * TotalArea / Sin(r(22.5)));
      Str (EquivSum:7:0, TextItem);
      OutTextXY (XCentre + Round(TotalRadius) + 10, YCentre, TextItem);
      SetFillStyle (EmptyFill, 0);
      Ellipse     (XCentre, YCentre, 0, 360,
                   Round(TotalRadius), Round(TotalRadius * Aspect));
      FillEllipse (XCentre, YCentre,
                   Round(TotalRadius), Round(TotalRadius * Aspect));
      EquivSum    := EquivSum * 10;
   end;
   ReadLn;
   CloseGraph;
end.