<template>
  <div id="quiz">
    <Header
      :numCorrect="numCorrect"
      :numTotal="numTotal"
      :field="field"
    />
    <b-container class="bv-example-row">
      <b-row>
        <b-col sm="12" offset="0.5">
          <QuestionBox
            v-if="questions.length"
            :currentQuestion="questions[index]"
            :next="next"
            :increment="increment"
          />
        </b-col>
      </b-row>
    </b-container>
  </div>
</template>

<script>
import Header from '../components/Header.vue'
import QuestionBox from '../components/QuestionBox.vue'

export default {
  name: 'Quiz',
  components: {
    Header,
    QuestionBox
  },
  data() {
    return {
      questions: [],
      index: 0,
      numCorrect: 0,
      numTotal: 0,
    }
  },
  props: [
    'field',
    'url'
  ],
  methods: {
    next() {
      this.index++
    },
    increment(isCorrect) {
      if (isCorrect) {
        this.numCorrect++
      }
      this.numTotal++
    }
  },
  mounted: function() {
    console.log(this.url)
    fetch(this.url, {
      method: 'get'
    })
    .then((response) => {
      return response.json()
    })
    .then((jsonData) => {
      this.questions = jsonData.results
    })
    console.log("Questions fetched")
  }
}
</script>

<style>
#quiz {
  font-family: Avenir, Helvetica, Arial, sans-serif;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  text-align: center;
  color: #2c3e50;
  margin-top: 60px;
}
</style>
