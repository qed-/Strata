/*
 * Copyright (C) 2021 - present by OpenGamma Inc. and the OpenGamma group of companies
 *
 * Please see distribution for license.
 */
package com.opengamma.strata.pricer.capfloor;

import static com.opengamma.strata.market.model.SabrParameterType.ALPHA;
import static com.opengamma.strata.market.model.SabrParameterType.BETA;
import static com.opengamma.strata.market.model.SabrParameterType.NU;
import static com.opengamma.strata.market.model.SabrParameterType.RHO;

import java.time.LocalDate;
import java.time.ZoneId;
import java.util.List;

import com.opengamma.strata.basics.currency.Currency;
import com.opengamma.strata.basics.currency.CurrencyAmount;
import com.opengamma.strata.basics.value.ValueDerivatives;
import com.opengamma.strata.collect.array.DoubleArray;
import com.opengamma.strata.market.sensitivity.PointSensitivityBuilder;
import com.opengamma.strata.pricer.ZeroRateSensitivity;
import com.opengamma.strata.pricer.impl.option.BlackFormulaRepository;
import com.opengamma.strata.pricer.impl.rate.ForwardOvernightCompoundedRateComputationFn;
import com.opengamma.strata.pricer.impl.volatility.smile.SabrFormulaData;
import com.opengamma.strata.pricer.impl.volatility.smile.SabrInarrearsVolatilityFunction;
import com.opengamma.strata.pricer.rate.RatesProvider;
import com.opengamma.strata.product.capfloor.OvernightInarrearsCapletFloorletPeriod;
import com.opengamma.strata.product.common.PutCall;
import com.opengamma.strata.product.rate.OvernightCompoundedRateComputation;

/**
 * Pricer for in-arrears caplets and floorlets (Asian style options) in the SABR with effective parameters approach.
 */
public class SabrOvernightInarrearsCapletFloorletPeriodPricer {

  /**
   * The default function for Asian in-arrears effective parameters.
   */
  private final SabrInarrearsVolatilityFunction sabrInarrearsFunction;

  /**
   * Default implementation.
   */
  public static final SabrOvernightInarrearsCapletFloorletPeriodPricer DEFAULT =
      new SabrOvernightInarrearsCapletFloorletPeriodPricer(SabrInarrearsVolatilityFunction.DEFAULT);

  /**
   * Creates an instance.
   * 
   * @param sabrInarrearsFunction  the function for Asian in-arrears effective parameters
   */
  private SabrOvernightInarrearsCapletFloorletPeriodPricer(SabrInarrearsVolatilityFunction sabrInarrearsFunction) {
    this.sabrInarrearsFunction = sabrInarrearsFunction;
  }
  
  /**
   * 
   * @param sabrInarrearsFunction
   * @return
   */
  public static SabrOvernightInarrearsCapletFloorletPeriodPricer of(SabrInarrearsVolatilityFunction sabrInarrearsFunction) {
    return new SabrOvernightInarrearsCapletFloorletPeriodPricer(sabrInarrearsFunction);
  }
  
  private static final ForwardOvernightCompoundedRateComputationFn ON_FUNCT = 
      ForwardOvernightCompoundedRateComputationFn.DEFAULT;

  /**
   * Computes the present value in the SABR model with effective parameters.
   * 
   * @param period  the caplet/floorlet period
   * @param ratesProvider  the rates provider
   * @param sabrVolatilities  the SABR volatilities
   * @return the present value
   */
  public CurrencyAmount presentValue(
      OvernightInarrearsCapletFloorletPeriod period,
      RatesProvider ratesProvider,
      SabrParametersIborCapletFloorletVolatilities sabrVolatilities) {

    Currency currency = period.getCurrency();
    if (ratesProvider.getValuationDate().isAfter(period.getPaymentDate())) {
      return CurrencyAmount.of(currency, 0d);
    }
    OvernightCompoundedRateComputation onComputation = period.getOvernightRate();
    LocalDate startDate = onComputation.getStartDate();
    LocalDate endDate = onComputation.getEndDate();
    double startTime = sabrVolatilities.relativeTime(startDate.atStartOfDay(ZoneId.of("Europe/Brussels"))); // xxx
    double endTime = sabrVolatilities.relativeTime(endDate.atStartOfDay(ZoneId.of("Europe/Brussels"))); // xxx
    double df = ratesProvider.discountFactor(currency, period.getPaymentDate());
    PutCall putCall = period.getPutCall();
    double strike = period.getStrike();
    double forward = ON_FUNCT
        .rate(onComputation, onComputation.getStartDate(), onComputation.getEndDate(), ratesProvider);
    double alpha = sabrVolatilities.alpha(startTime); // parameters at start of composition period, for coherence with term rate caplets
    double beta = sabrVolatilities.beta(startTime);
    double rho = sabrVolatilities.rho(startTime);
    double nu = sabrVolatilities.nu(startTime);
    double shift = sabrVolatilities.shift(startTime);
    SabrFormulaData sabr = SabrFormulaData.of(alpha, beta, rho, nu);
    SabrFormulaData sabrAdjusted = sabrInarrearsFunction.effectiveSabr(sabr, startTime, endTime);
    double volatility = sabrVolatilities.getParameters().getSabrVolatilityFormula()
        .volatility(forward + shift, strike + shift, endTime,
            sabrAdjusted.getAlpha(), sabrAdjusted.getBeta(), sabrAdjusted.getRho(), sabrAdjusted.getNu());
    double price = df * period.getYearFraction() *
        BlackFormulaRepository.price(forward + shift, strike + shift, endTime, volatility, putCall.isCall());
    return CurrencyAmount.of(currency, price * period.getNotional());
  }

  /**
   * Computes the present value sensitivity to the rate with "sticky SABR model parameters".
   * 
   * @param period  the caplet/floorlet period
   * @param ratesProvider  the rates provider
   * @param sabrVolatilities  the SABR volatilities
   * @return the present value rate sensitivity
   */
  public PointSensitivityBuilder presentValueSensitivityRatesStickyModel(
      OvernightInarrearsCapletFloorletPeriod period,
      RatesProvider ratesProvider,
      SabrParametersIborCapletFloorletVolatilities volatilities) {

    Currency currency = period.getCurrency();
    if (ratesProvider.getValuationDate().isAfter(period.getPaymentDate())) {
      return PointSensitivityBuilder.none();
    }
    OvernightCompoundedRateComputation onComputation = period.getOvernightRate();
    LocalDate startDate = onComputation.getStartDate();
    LocalDate endDate = onComputation.getEndDate();
    double startTime = volatilities.relativeTime(startDate.atStartOfDay(ZoneId.of("Europe/Brussels"))); // xxx
    double endTime = volatilities.relativeTime(endDate.atStartOfDay(ZoneId.of("Europe/Brussels"))); // xxx
    double df = ratesProvider.discountFactor(currency, period.getPaymentDate());
    PutCall putCall = period.getPutCall();
    double strike = period.getStrike();
    double forward = ON_FUNCT
        .rate(onComputation, onComputation.getStartDate(), onComputation.getEndDate(), ratesProvider);
    double alpha = volatilities.alpha(startTime);
    double beta = volatilities.beta(startTime);
    double rho = volatilities.rho(startTime);
    double nu = volatilities.nu(startTime);
    double shift = volatilities.shift(startTime);
    SabrFormulaData sabr = SabrFormulaData.of(alpha, beta, rho, nu);
    SabrFormulaData sabrAdjusted = sabrInarrearsFunction.effectiveSabr(sabr, startTime, endTime);
    ValueDerivatives volatility = volatilities.getParameters().getSabrVolatilityFormula()
        .volatilityAdjoint(forward + shift, strike + shift, endTime,
            sabrAdjusted.getAlpha(), sabrAdjusted.getBeta(), sabrAdjusted.getRho(), sabrAdjusted.getNu());
    ValueDerivatives price = BlackFormulaRepository
        .priceAdjoint(forward + shift, strike + shift, endTime, volatility.getValue(), putCall.isCall());
    double pv = price.getValue() * df * period.getYearFraction() * period.getNotional();
    // Backward sweep
    double pvBar = 1.0;
    double priceBar =  df * period.getYearFraction() * period.getNotional() * pvBar;
    double dfBar = pv / df * pvBar;
    double forwardBar = price.getDerivative(0) * priceBar;
    double volatilityBar = price.getDerivative(3) * priceBar;
    forwardBar += volatility.getDerivative(0) * volatilityBar;
    PointSensitivityBuilder dforwarddr = ON_FUNCT
        .rateSensitivity(onComputation, onComputation.getStartDate(), onComputation.getEndDate(), ratesProvider);
    ZeroRateSensitivity ddfdr = ratesProvider
        .discountFactors(currency).zeroRatePointSensitivity(period.getPaymentDate());
    return ddfdr.multipliedBy(dfBar).combinedWith(dforwarddr.multipliedBy(forwardBar));
  }

  /**
   * Computes the present value sensitivity to the SABR model parameters.
   * 
   * @param period  the caplet/floorlet period
   * @param ratesProvider  the rates provider
   * @param sabrVolatilities  the SABR volatilities
   * @return the present value model parameters sensitivity
   */
  public PointSensitivityBuilder presentValueSensitivityModelParamsSabr(
      OvernightInarrearsCapletFloorletPeriod period,
      RatesProvider ratesProvider,
      SabrParametersIborCapletFloorletVolatilities volatilities) {

    Currency currency = period.getCurrency();
    if (ratesProvider.getValuationDate().isAfter(period.getPaymentDate())) {
      return PointSensitivityBuilder.none();
    }
    OvernightCompoundedRateComputation onComputation = period.getOvernightRate();
    LocalDate startDate = onComputation.getStartDate();
    LocalDate endDate = onComputation.getEndDate();
    double startTime = volatilities.relativeTime(startDate.atStartOfDay(ZoneId.of("Europe/Brussels"))); // xxx
    double endTime = volatilities.relativeTime(endDate.atStartOfDay(ZoneId.of("Europe/Brussels"))); // xxx
    double df = ratesProvider.discountFactor(currency, period.getPaymentDate());
    PutCall putCall = period.getPutCall();
    double strike = period.getStrike();
    double forward = ON_FUNCT
        .rate(onComputation, onComputation.getStartDate(), onComputation.getEndDate(), ratesProvider);
    double alpha = volatilities.alpha(startTime);
    double beta = volatilities.beta(startTime);
    double rho = volatilities.rho(startTime);
    double nu = volatilities.nu(startTime);
    double shift = volatilities.shift(startTime);
    SabrFormulaData sabr = SabrFormulaData.of(alpha, beta, rho, nu);
    List<ValueDerivatives> sabrAdjusted = sabrInarrearsFunction.effectiveSabrAd(sabr, startTime, endTime);
    ValueDerivatives volatility = volatilities.getParameters().getSabrVolatilityFormula()
        .volatilityAdjoint(forward + shift, strike + shift, endTime,
            sabrAdjusted.get(0).getValue(), sabrAdjusted.get(1).getValue(),
            sabrAdjusted.get(2).getValue(), sabrAdjusted.get(3).getValue());
    ValueDerivatives price = BlackFormulaRepository
        .priceAdjoint(forward + shift, strike + shift, endTime, volatility.getValue(), putCall.isCall());
//    double pv = price.getValue() * df * period.getYearFraction() * period.getNotional();
    // Backward sweep
    double pvBar = 1.0;
    double priceBar = df * period.getYearFraction() * period.getNotional() * pvBar;
    double volatilityBar = price.getDerivative(3) * priceBar;
    double alphaHatBar = volatility.getDerivative(2) * volatilityBar;
    double betaHatBar = volatility.getDerivative(3) * volatilityBar;
    double rhoHatBar = volatility.getDerivative(4) * volatilityBar;
    double nuHatBar = volatility.getDerivative(5) * volatilityBar;
    DoubleArray paramHat = sabrAdjusted.get(0).getDerivatives().multipliedBy(alphaHatBar);
    paramHat = paramHat.plus(sabrAdjusted.get(1).getDerivatives().multipliedBy(betaHatBar));
    paramHat = paramHat.plus(sabrAdjusted.get(2).getDerivatives().multipliedBy(rhoHatBar));
    paramHat = paramHat.plus(sabrAdjusted.get(3).getDerivatives().multipliedBy(nuHatBar));

    IborCapletFloorletVolatilitiesName name = volatilities.getName();
    return PointSensitivityBuilder.of(
        IborCapletFloorletSabrSensitivity.of(name, startTime, ALPHA, currency, paramHat.get(0)),
        IborCapletFloorletSabrSensitivity.of(name, startTime, BETA, currency, paramHat.get(1)),
        IborCapletFloorletSabrSensitivity.of(name, startTime, RHO, currency, paramHat.get(2)),
        IborCapletFloorletSabrSensitivity.of(name, startTime, NU, currency, paramHat.get(3)));
  }

}
