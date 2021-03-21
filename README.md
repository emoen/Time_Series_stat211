# Time_Series_stat211
This course gives an introduction to linear time series models, such as autoregressive, moving average and ARMA models. Moreover, it is shown how the empirical autocorrelation and partial correlation can be used to identify the model. The Durbin- Levinson, the innovation algorithm and the theory for optimal forecasts are explained. The last part of the course gives an introduction to methods of estimation. Empirical modelling using the AIC and FPE criteria is mentioned as is ARCH and GARCH models.

## ARMA(1) model on deaths dataset
![deaths_detrended_deseasoned_w_prediction.png](https://github.com/emoen/Time_Series_stat211/blob/main/hw1/deaths_detrended_deseasoned_w_prediction.png)

## Simulating AR(1) model



AR(1) with ø=0.7, plotting Auto Correlation(ACF) and Partial Auto Correlation Function(PACF). The ACF tails off while the PACF cuts off after one lag
![ACF_PACF_AR_1_b.png](https://github.com/emoen/Time_Series_stat211/blob/main/hw3/ACF_PACF_AR_1_b.png)


|f\Mdl|      AR(P)     | MA(q)          | ARMA(p,q)|
|-----|----------------|----------------|----------|
|ACF   |Tails off      | cuts off at=q  | Tails off|
|PACF  |cuts off at t=p| Tails off      | Tails off|

## Simulating AR(2) model
Draw the rectangle defined by {φ2 = 1, φ2 − φ1 = 1} in a φ1φ2-coordinate system:

![AR_2_stationarity.png](https://github.com/emoen/Time_Series_stat211/blob/main/hw4/AR_2_stationarity.png)


