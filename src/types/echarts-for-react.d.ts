// echarts-for-react@3.x was built against React 15-18 class component types.
// React 19 changed the Component base, breaking JSX usage. This patches the
// exported class to be compatible with React 19's JSX element type check.
import type { EChartsReactProps } from 'echarts-for-react';
import type { Component } from 'react';

declare module 'echarts-for-react' {
  class EChartsReact extends Component<EChartsReactProps> {
    getEchartsInstance(): import('echarts').ECharts;
  }
  export default EChartsReact;
}
