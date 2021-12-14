/*
 * Copyright (C) 2021 - present by OpenGamma Inc. and the OpenGamma group of companies
 *
 * Please see distribution for license.
 */
package com.opengamma.strata.pricer.capfloor;

import static com.opengamma.strata.basics.date.HolidayCalendarIds.EUTA;
import static com.opengamma.strata.basics.currency.Currency.EUR;
import static com.opengamma.strata.basics.index.IborIndices.EUR_EURIBOR_3M; // TODO: replace with a ON benchmark
import static com.opengamma.strata.basics.index.OvernightIndices.EUR_ESTR;
import static com.opengamma.strata.collect.TestHelper.dateUtc;
import static org.assertj.core.api.Assertions.assertThat;

import java.time.LocalDate;
import java.time.ZoneId;
import java.time.ZonedDateTime;

import org.assertj.core.data.Offset;
import org.junit.jupiter.api.Test;

import com.opengamma.strata.basics.ReferenceData;
import com.opengamma.strata.basics.currency.CurrencyAmount;
import com.opengamma.strata.basics.date.DayCounts;
import com.opengamma.strata.basics.date.HolidayCalendar;
import com.opengamma.strata.basics.index.OvernightIndexObservation;
import com.opengamma.strata.collect.timeseries.LocalDateDoubleTimeSeries;
import com.opengamma.strata.market.curve.CurveName;
import com.opengamma.strata.market.param.CurrencyParameterSensitivities;
import com.opengamma.strata.market.param.CurrencyParameterSensitivity;
import com.opengamma.strata.market.sensitivity.PointSensitivities;
import com.opengamma.strata.pricer.ZeroRateDiscountFactors;
import com.opengamma.strata.pricer.impl.option.BlackFormulaRepository;
import com.opengamma.strata.pricer.impl.option.NormalFormulaRepository;
import com.opengamma.strata.pricer.impl.swap.DiscountingRatePaymentPeriodPricer;
import com.opengamma.strata.pricer.impl.volatility.smile.SabrFormulaData;
import com.opengamma.strata.pricer.impl.volatility.smile.SabrHaganVolatilityFunctionProvider;
import com.opengamma.strata.pricer.impl.volatility.smile.SabrInarrearsVolatilityFunction;
import com.opengamma.strata.pricer.rate.DiscountOvernightIndexRates;
import com.opengamma.strata.pricer.rate.ImmutableRatesProvider;
import com.opengamma.strata.pricer.sensitivity.RatesFiniteDifferenceSensitivityCalculator;
import com.opengamma.strata.product.capfloor.IborCapletFloorletPeriod;
import com.opengamma.strata.product.capfloor.OvernightInarrearsCapletFloorletPeriod;
import com.opengamma.strata.product.common.PutCall;
import com.opengamma.strata.product.rate.FixedRateComputation;
import com.opengamma.strata.product.rate.IborRateComputation;
import com.opengamma.strata.product.rate.OvernightCompoundedRateComputation;
import com.opengamma.strata.product.swap.RateAccrualPeriod;
import com.opengamma.strata.product.swap.RatePaymentPeriod;

/**
 * Test {@link SabrOvernightInarrearsCapletFloorletPeriodPricer}.
 */
public class SabrOvernightInarrearsCapletFloorletPeriodPricerTest {

  private static final ReferenceData REF_DATA = ReferenceData.standard();
  private static final HolidayCalendar EUTA_IMPL = REF_DATA.getValue(EUTA);
  private static final ZonedDateTime VALUATION = dateUtc(2021, 12, 20);
  private static final LocalDate START_DATE = LocalDate.of(2022, 6, 22);
  private static final LocalDate END_DATE = LocalDate.of(2022, 9, 22);
  private static final LocalDate PAYMENT_DATE = END_DATE.plusMonths(1);
  private static final double NOTIONAL = 1_000_000.0d;
  private static final double STRIKE = 0.0155;
  private static final double ACCRUAL_FACTOR = 0.30;
  private static final OvernightCompoundedRateComputation RATE_COMP =
      OvernightCompoundedRateComputation.of(EUR_ESTR, START_DATE, END_DATE, REF_DATA);

  private static final OvernightInarrearsCapletFloorletPeriod CAPLET_LONG =
      OvernightInarrearsCapletFloorletPeriod.builder()
          .caplet(STRIKE)
          .startDate(START_DATE)
          .endDate(END_DATE)
          .paymentDate(PAYMENT_DATE)
          .yearFraction(ACCRUAL_FACTOR)
          .notional(NOTIONAL)
          .overnightRate(RATE_COMP)
          .build();
  private static final OvernightInarrearsCapletFloorletPeriod CAPLET_SHORT =
      OvernightInarrearsCapletFloorletPeriod.builder()
          .caplet(STRIKE)
          .startDate(START_DATE)
          .endDate(END_DATE)
          .paymentDate(PAYMENT_DATE)
          .yearFraction(ACCRUAL_FACTOR)
          .notional(-NOTIONAL)
          .overnightRate(RATE_COMP)
          .build();
  private static final OvernightInarrearsCapletFloorletPeriod FLOORLET_LONG =
      OvernightInarrearsCapletFloorletPeriod.builder()
          .floorlet(STRIKE)
          .startDate(START_DATE)
          .endDate(END_DATE)
          .paymentDate(PAYMENT_DATE)
          .yearFraction(ACCRUAL_FACTOR)
          .notional(NOTIONAL)
          .overnightRate(RATE_COMP)
          .build();
  private static final OvernightInarrearsCapletFloorletPeriod FLOORLET_SHORT =
      OvernightInarrearsCapletFloorletPeriod.builder()
          .floorlet(STRIKE)
          .startDate(START_DATE)
          .endDate(END_DATE)
          .paymentDate(PAYMENT_DATE)
          .yearFraction(ACCRUAL_FACTOR)
          .notional(-NOTIONAL)
          .overnightRate(RATE_COMP)
          .build();

  // valuation date before fixing date
  private static final ImmutableRatesProvider RATES = IborCapletFloorletSabrRateVolatilityDataSet.getRatesProvider(
      VALUATION.toLocalDate(), EUR_EURIBOR_3M, LocalDateDoubleTimeSeries.empty());
  private static final SabrParametersIborCapletFloorletVolatilities VOLS = IborCapletFloorletSabrRateVolatilityDataSet
      .getVolatilities(VALUATION, EUR_EURIBOR_3M);

  private static final SabrIborCapletFloorletPeriodPricer PRICER_IBOR = SabrIborCapletFloorletPeriodPricer.DEFAULT;
  private static final SabrOvernightInarrearsCapletFloorletPeriodPricer PRICER_ON_INARREARS_Q1 =
      SabrOvernightInarrearsCapletFloorletPeriodPricer.DEFAULT;
//  private static final SabrInarrearsVolatilityFunction FUNCTION_DEFAULT =
//      SabrInarrearsVolatilityFunction.DEFAULT;
  private static final double Q_OTHER = 1.1;
  private static final SabrInarrearsVolatilityFunction FUNCTION_OTHER =
      SabrInarrearsVolatilityFunction.of(Q_OTHER);
  private static final SabrOvernightInarrearsCapletFloorletPeriodPricer PRICER_ON_INARREARS_QOTHER =
      SabrOvernightInarrearsCapletFloorletPeriodPricer.of(FUNCTION_OTHER);
  private static final SabrHaganVolatilityFunctionProvider SABR_FORMULA =
      SabrHaganVolatilityFunctionProvider.DEFAULT;
  private static final double EPS_FD = 1.0e-6;
  private static final RatesFiniteDifferenceSensitivityCalculator FD_CAL =
      new RatesFiniteDifferenceSensitivityCalculator(EPS_FD);
  private static final DiscountingRatePaymentPeriodPricer PRICER_PERIOD =
      DiscountingRatePaymentPeriodPricer.DEFAULT;

  private static final Offset<Double> TOLERANCE_SMALL_IV = Offset.offset(1.0E-4);
  private static final Offset<Double> TOLERANCE_SMALL_PV = Offset.offset(1.0E+1);
  private static final Offset<Double> TOLERANCE_PV = Offset.offset(1.0E-4);
  private static final Offset<Double> TOLERANCE_VEGA = Offset.offset(1.0E-2);
  private static final Double TOLERANCE_DELTA = 2.0E+1; // 0.002 / bps on 1m

  /* The present value of the in-arrears option is higher than the in-advance. */
  @Test
  public void presentValue_higher() {
    ZeroRateDiscountFactors onDf =
        (ZeroRateDiscountFactors) ((DiscountOvernightIndexRates) RATES.overnightIndexRates(EUR_ESTR))
            .getDiscountFactors();
    ImmutableRatesProvider ratesOn = RATES.toBuilder()
        .iborIndexCurve(EUR_EURIBOR_3M, onDf.getCurve()).build(); // Change IBOR curve to ON for comparison
    IborRateComputation iborComp = IborRateComputation.of(EUR_EURIBOR_3M, LocalDate.of(2022, 6, 20), REF_DATA);
    IborCapletFloorletPeriod iborCap = IborCapletFloorletPeriod.builder()
        .caplet(STRIKE)
        .startDate(START_DATE)
        .endDate(END_DATE)
        .paymentDate(PAYMENT_DATE)
        .yearFraction(ACCRUAL_FACTOR)
        .notional(NOTIONAL)
        .iborRate(iborComp)
        .build();
    CurrencyAmount inAdvancePv = PRICER_IBOR.presentValue(iborCap, ratesOn, VOLS);
    CurrencyAmount inArrearsPv = PRICER_ON_INARREARS_Q1.presentValue(CAPLET_LONG, ratesOn, VOLS);
    assertThat(inAdvancePv.getAmount() < inArrearsPv.getAmount()).isTrue();
  }

  /* The implied volatility of the in-arrears option converges to the in-advance when the maturity tends to infinity. */
  @Test
  public void impliedvol_converge() {
    ZeroRateDiscountFactors onDf =
        (ZeroRateDiscountFactors) ((DiscountOvernightIndexRates) RATES.overnightIndexRates(EUR_ESTR))
            .getDiscountFactors();
    ImmutableRatesProvider ratesOn = RATES.toBuilder()
        .iborIndexCurve(EUR_EURIBOR_3M, onDf.getCurve()).build(); // Change IBOR curve to ON for comparison
    for (int i = 50; i <= 100; i += 10) {
      LocalDate fixingDate = EUTA_IMPL.nextOrSame(START_DATE.plusYears(i));
      IborRateComputation iborComp = IborRateComputation.of(EUR_EURIBOR_3M, fixingDate, REF_DATA);
      LocalDate startDate = iborComp.getEffectiveDate();
      LocalDate endDate = iborComp.getMaturityDate();
      OvernightCompoundedRateComputation onComp =
          OvernightCompoundedRateComputation.of(EUR_ESTR, startDate, endDate, REF_DATA);
      OvernightIndexObservation onObs =
          OvernightIndexObservation.of(EUR_ESTR, startDate, REF_DATA);
      OvernightInarrearsCapletFloorletPeriod inArrearsCap = OvernightInarrearsCapletFloorletPeriod.builder()
          .caplet(STRIKE)
          .startDate(startDate)
          .endDate(endDate)
          .yearFraction(ACCRUAL_FACTOR)
          .notional(NOTIONAL)
          .overnightRate(onComp)
          .build();
      IborCapletFloorletPeriod iborCap = IborCapletFloorletPeriod.builder()
          .caplet(STRIKE)
          .startDate(startDate)
          .endDate(endDate)
          .yearFraction(ACCRUAL_FACTOR)
          .notional(NOTIONAL)
          .iborRate(iborComp)
          .build();
      double forward = ratesOn.overnightIndexRates(EUR_ESTR).periodRate(onObs, endDate);
      double num = NOTIONAL * ACCRUAL_FACTOR * onDf.discountFactor(endDate);
      double timeToExpiry = VOLS.relativeTime(fixingDate.atStartOfDay(ZoneId.of("Europe/Brussels")));
      CurrencyAmount inArrearsPv = PRICER_ON_INARREARS_Q1.presentValue(inArrearsCap, ratesOn, VOLS);
      double inArrearsIv = NormalFormulaRepository
          .impliedVolatility(inArrearsPv.getAmount(), forward, STRIKE, timeToExpiry, 0.01, num, PutCall.CALL);
      CurrencyAmount inAdvancePv = PRICER_IBOR.presentValue(iborCap, ratesOn, VOLS);
      double inAdvanceIv = NormalFormulaRepository
          .impliedVolatility(inAdvancePv.getAmount(), forward, STRIKE, timeToExpiry, 0.01, num, PutCall.CALL);
      assertThat(inArrearsIv).isEqualTo(inAdvanceIv, TOLERANCE_SMALL_IV);
//      System.out.println("Difference (" + i + " years): " + (inArrearsIv - inAdvanceIv));
      // TODO: remove
    }
  }

  /* The value of the in-arrears option converges to the in-advance when q tends to infinity. */
  @Test
  public void presentvalue_converge_pinfinity() {
    ZeroRateDiscountFactors onDf =
        (ZeroRateDiscountFactors) ((DiscountOvernightIndexRates) RATES.overnightIndexRates(EUR_ESTR))
            .getDiscountFactors();
    ImmutableRatesProvider ratesOn = RATES.toBuilder()
        .iborIndexCurve(EUR_EURIBOR_3M, onDf.getCurve()).build(); // Change IBOR curve to ON for comparison
    LocalDate fixingDate = LocalDate.of(2022, 6, 20);
    IborRateComputation iborComp = IborRateComputation.of(EUR_EURIBOR_3M,fixingDate , REF_DATA);
    IborCapletFloorletPeriod iborCap = IborCapletFloorletPeriod.builder()
        .caplet(STRIKE)
        .startDate(START_DATE)
        .endDate(END_DATE)
        .paymentDate(PAYMENT_DATE)
        .yearFraction(ACCRUAL_FACTOR)
        .notional(NOTIONAL)
        .iborRate(iborComp)
        .build();
    OvernightIndexObservation onObs =
        OvernightIndexObservation.of(EUR_ESTR, START_DATE, REF_DATA);
    double forward = ratesOn.overnightIndexRates(EUR_ESTR).periodRate(onObs, END_DATE);
    double num = NOTIONAL * ACCRUAL_FACTOR * onDf.discountFactor(END_DATE);
    double timeToExpiry = VOLS.relativeTime(fixingDate.atStartOfDay(ZoneId.of("Europe/Brussels")));
    for (int i = 0; i < 10; i++) {
      double q = 100 + 10 * i; // Exaggerated q value to see the convergence
      SabrInarrearsVolatilityFunction functionQ = SabrInarrearsVolatilityFunction.of(q);
      SabrOvernightInarrearsCapletFloorletPeriodPricer pricerQ =
          SabrOvernightInarrearsCapletFloorletPeriodPricer.of(functionQ);
      double inAdvancePv = PRICER_IBOR.presentValue(iborCap, ratesOn, VOLS).getAmount();
      double inArrearsPv = pricerQ.presentValue(CAPLET_LONG, ratesOn, VOLS).getAmount();
      double inAdvanceIv = NormalFormulaRepository
          .impliedVolatility(inAdvancePv, forward, STRIKE, timeToExpiry, 0.01, num, PutCall.CALL);
      double inArrearsIv = NormalFormulaRepository
          .impliedVolatility(inArrearsPv, forward, STRIKE, timeToExpiry, 0.01, num, PutCall.CALL);
//    System.out.println("Difference (" + i + " years): " + (inArrearsIv - inAdvanceIv));
      assertThat(inArrearsIv).isEqualTo(inAdvanceIv, TOLERANCE_SMALL_IV);
    }
  }

  @Test
  public void presentValue_before_parity() {
    CurrencyAmount pvCapLongComputed = PRICER_ON_INARREARS_QOTHER.presentValue(CAPLET_LONG, RATES, VOLS);
    CurrencyAmount pvFloorShortComputed = PRICER_ON_INARREARS_QOTHER.presentValue(FLOORLET_SHORT, RATES, VOLS);
    CurrencyAmount pvCapShortComputed = PRICER_ON_INARREARS_QOTHER.presentValue(CAPLET_SHORT, RATES, VOLS);
    CurrencyAmount pvFloorLongComputed = PRICER_ON_INARREARS_QOTHER.presentValue(FLOORLET_LONG, RATES, VOLS);
    RateAccrualPeriod onAccrual = RateAccrualPeriod.builder()
        .rateComputation(RATE_COMP)
        .startDate(START_DATE)
        .endDate(END_DATE)
        .yearFraction(ACCRUAL_FACTOR).build();
    RatePaymentPeriod onPeriod = RatePaymentPeriod.builder()
        .paymentDate(PAYMENT_DATE)
        .accrualPeriods(onAccrual)
        .dayCount(DayCounts.ACT_360)
        .currency(EUR)
        .notional(NOTIONAL).build();
    RateAccrualPeriod fixedAccrual = RateAccrualPeriod.builder()
        .rateComputation(FixedRateComputation.of(STRIKE))
        .startDate(START_DATE)
        .endDate(END_DATE)
        .yearFraction(ACCRUAL_FACTOR).build();
    RatePaymentPeriod fixedPeriod = RatePaymentPeriod.builder()
        .paymentDate(PAYMENT_DATE)
        .accrualPeriods(fixedAccrual)
        .dayCount(DayCounts.ACT_360)
        .currency(EUR)
        .notional(NOTIONAL).build();
    double pvOnPeriod = PRICER_PERIOD.presentValue(onPeriod, RATES);
    double pvFixedPeriod = PRICER_PERIOD.presentValue(fixedPeriod, RATES);
    assertThat(pvCapLongComputed.getAmount() + pvFloorShortComputed.getAmount() + pvFixedPeriod)
        .isEqualTo(pvOnPeriod, TOLERANCE_SMALL_IV);
    assertThat(pvCapShortComputed.getAmount() + pvFloorLongComputed.getAmount() - pvFixedPeriod)
        .isEqualTo(-pvOnPeriod, TOLERANCE_SMALL_IV);
  }

  @Test
  public void presentValue_before_formula() {
    CurrencyAmount pvLongComputed = PRICER_ON_INARREARS_QOTHER.presentValue(CAPLET_LONG, RATES, VOLS);
    OvernightIndexObservation onObs = OvernightIndexObservation.of(EUR_ESTR, START_DATE, REF_DATA);
    double startTime = VOLS.relativeTime(START_DATE.atStartOfDay(ZoneId.of("Europe/Brussels")));
    double endTime = VOLS.relativeTime(END_DATE.atStartOfDay(ZoneId.of("Europe/Brussels")));
    double alpha = VOLS.alpha(startTime);
    double beta = VOLS.beta(startTime);
    double rho = VOLS.rho(startTime);
    double nu = VOLS.nu(startTime);
    double shift = VOLS.shift(startTime);
    SabrFormulaData sabr = SabrFormulaData.of(alpha, beta, rho, nu);
    SabrFormulaData sabrEffective = FUNCTION_OTHER.effectiveSabr(sabr, startTime, endTime);
    double forward = RATES.overnightIndexRates(EUR_ESTR).periodRate(onObs, END_DATE);
    double dfPayment = RATES.discountFactor(EUR, PAYMENT_DATE);
    double volatility = SABR_FORMULA.volatility(forward + shift, STRIKE + shift, endTime, sabrEffective);
    double pvExpected = dfPayment * ACCRUAL_FACTOR * NOTIONAL *
        BlackFormulaRepository.price(forward + shift, STRIKE + shift, endTime, volatility, true);
    assertThat(pvExpected).isEqualTo(pvLongComputed.getAmount(), TOLERANCE_PV);
    CurrencyAmount pvShortComputed = PRICER_ON_INARREARS_QOTHER.presentValue(CAPLET_SHORT, RATES, VOLS);
    assertThat(pvShortComputed.getAmount()).isEqualTo(-pvLongComputed.getAmount(), TOLERANCE_PV);
  }

  @Test
  public void presentValue_before_rate_sensitivity() {
    PointSensitivities ptsCapletLongComputed = PRICER_ON_INARREARS_QOTHER
        .presentValueSensitivityRatesStickyModel(CAPLET_LONG, RATES, VOLS).build();
    CurrencyParameterSensitivities psCapletLongComputed = RATES.parameterSensitivity(ptsCapletLongComputed);
    CurrencyParameterSensitivities psCapletLongExpected =
        FD_CAL.sensitivity(RATES, p -> PRICER_ON_INARREARS_QOTHER.presentValue(CAPLET_LONG, p, VOLS));
    assertThat(psCapletLongComputed.equalWithTolerance(psCapletLongExpected, TOLERANCE_DELTA)).isTrue();
    PointSensitivities ptsFloorletShortComputed = PRICER_ON_INARREARS_QOTHER
        .presentValueSensitivityRatesStickyModel(FLOORLET_SHORT, RATES, VOLS).build();
    CurrencyParameterSensitivities psFloorletShortComputed = RATES.parameterSensitivity(ptsFloorletShortComputed);
    CurrencyParameterSensitivities psFloorletShortExpected =
        FD_CAL.sensitivity(RATES, p -> PRICER_ON_INARREARS_QOTHER.presentValue(FLOORLET_SHORT, p, VOLS));
    assertThat(psFloorletShortComputed.equalWithTolerance(psFloorletShortExpected, TOLERANCE_DELTA)).isTrue();
  }

  @Test
  public void presentValue_before_param_sensitivity() {
    PointSensitivities ptsCapletLongComputed = PRICER_ON_INARREARS_QOTHER
        .presentValueSensitivityModelParamsSabr(CAPLET_LONG, RATES, VOLS).build();
    CurrencyParameterSensitivities ps = VOLS.parameterSensitivity(ptsCapletLongComputed);
    CurrencyParameterSensitivity psAlpha = ps.getSensitivity(CurveName.of("Test-SABR-Alpha"), EUR);
    CurrencyParameterSensitivity psBeta = ps.getSensitivity(CurveName.of("Test-SABR-Beta"), EUR);
    CurrencyParameterSensitivity psRho = ps.getSensitivity(CurveName.of("Test-SABR-Rho"), EUR);
    CurrencyParameterSensitivity psNu = ps.getSensitivity(CurveName.of("Test-SABR-Nu"), EUR);
    double shift = 1.0E-6;
    double pvStart = PRICER_ON_INARREARS_QOTHER.presentValue(CAPLET_LONG, RATES, VOLS).getAmount();
    for (int i = 0; i < 6; i++) { // size of Alpha curve
      SabrParametersIborCapletFloorletVolatilities volsAlphaShifted = IborCapletFloorletSabrRateVolatilityDataSet
          .getVolatilitiesShiftAlpha(VALUATION, EUR_EURIBOR_3M, i, shift);
      double pvAlphaShifted = PRICER_ON_INARREARS_QOTHER.presentValue(CAPLET_LONG, RATES, volsAlphaShifted).getAmount();
      double psAlphaExpected = (pvAlphaShifted - pvStart) / shift;
      assertThat(psAlphaExpected).isEqualTo(psAlpha.getSensitivity().get(i), TOLERANCE_VEGA);
      SabrParametersIborCapletFloorletVolatilities volsBetaShifted = IborCapletFloorletSabrRateVolatilityDataSet
          .getVolatilitiesShiftBeta(VALUATION, EUR_EURIBOR_3M, i, shift);
      double pvBetaShifted = PRICER_ON_INARREARS_QOTHER.presentValue(CAPLET_LONG, RATES, volsBetaShifted).getAmount();
      double psBetaExpected = (pvBetaShifted - pvStart) / shift;
      assertThat(psBetaExpected).isEqualTo(psBeta.getSensitivity().get(i), TOLERANCE_VEGA);
      SabrParametersIborCapletFloorletVolatilities volsRhoShifted = IborCapletFloorletSabrRateVolatilityDataSet
          .getVolatilitiesShiftRho(VALUATION, EUR_EURIBOR_3M, i, shift);
      double pvRhoShifted = PRICER_ON_INARREARS_QOTHER.presentValue(CAPLET_LONG, RATES, volsRhoShifted).getAmount();
      double psRhoExpected = (pvRhoShifted - pvStart) / shift;
      assertThat(psRhoExpected).isEqualTo(psRho.getSensitivity().get(i), TOLERANCE_VEGA);
      SabrParametersIborCapletFloorletVolatilities volsNuShifted = IborCapletFloorletSabrRateVolatilityDataSet
          .getVolatilitiesShiftNu(VALUATION, EUR_EURIBOR_3M, i, shift);
      double pvNuShifted = PRICER_ON_INARREARS_QOTHER.presentValue(CAPLET_LONG, RATES, volsNuShifted).getAmount();
      double psNuExpected = (pvNuShifted - pvStart) / shift;
      assertThat(psNuExpected).isEqualTo(psNu.getSensitivity().get(i), TOLERANCE_VEGA);
    }
  }

  @Test
  public void presentValue_after_formula() {

  }

  // TODO: after start

}
