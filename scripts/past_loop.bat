for /l %%p in (1,1,25) do (
	echo Simulacion pasado: %%p
	gams NoSlack.gms --past=%%p
)
pause